
# ============================================================
# TRANSCRIPTOME-X
# STEP 8 - FUNCTIONAL ENRICHMENT (STRICT MAPPING)
# ============================================================
#
# Dataset: GSE266566
# Cell line: MDA-MB-231
# Biological context: NRP1 knockdown
#
# Mapping strategy:
# Strict one-to-one Ensembl-Entrez associations only.
#
# Analyses:
#   - GO Biological Process
#   - GO Molecular Function
#   - GO Cellular Component
#   - Reactome Pathways
#
# Separate analyses:
#   - UP-regulated genes
#   - DOWN-regulated genes
#   - siNRP1_3 vs siNT
#   - siNRP1_5 vs siNT
#
# Statistical threshold: FDR < 0.05
# Background: all edgeR-tested genes with strict mapping
#
# ============================================================


# 1. CONFIGURATION
# ============================================================

project_dir <- getwd()

input_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "differential_expression",
  "annotated"
)

previous_mapping_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "functional_enrichment",
  "mapping"
)

output_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "functional_enrichment_strict"
)

output_subdirs <- c(
  "results",
  "plots",
  "mapping",
  "reports"
)

for (subdir in output_subdirs) {
  dir.create(
    file.path(output_dir, subdir),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# 2. LOAD PACKAGES
# ============================================================

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "enrichplot",
  "AnnotationDbi",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(AnnotationDbi)
  library(ggplot2)
})


# 3. INPUT FILES
# ============================================================

input_files <- c(
  siNRP1_3_vs_siNT = "siNRP1_3_vs_siNT_annotated.csv",
  siNRP1_5_vs_siNT = "siNRP1_5_vs_siNT_annotated.csv"
)

de_data <- list()

for (comparison in names(input_files)) {
  
  input_path <- file.path(
    input_dir,
    input_files[[comparison]]
  )
  
  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path)
  }
  
  current_data <- read.csv(
    input_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  if (!all(required_columns %in% names(current_data))) {
    stop(
      "Missing required columns in ",
      input_files[[comparison]]
    )
  }
  
  if (anyDuplicated(current_data$gene_id)) {
    stop(
      "Duplicated gene IDs in ",
      input_files[[comparison]]
    )
  }
  
  if (anyNA(current_data$gene_id)) {
    stop(
      "Missing gene IDs in ",
      input_files[[comparison]]
    )
  }
  
  de_data[[comparison]] <- current_data
  
  cat(
    "\nComparison:",
    comparison,
    "| Genes:",
    nrow(current_data),
    "\n"
  )
}


# 4. VERIFY COMMON GENE UNIVERSE
# ============================================================

background_original <- sub(
  "\\.[0-9]+$",
  "",
  de_data[[1]]$gene_id
)

background_second <- sub(
  "\\.[0-9]+$",
  "",
  de_data[[2]]$gene_id
)

if (!setequal(background_original, background_second)) {
  stop(
    "The two comparisons do not share the same gene universe."
  )
}

if (anyDuplicated(background_original)) {
  stop(
    "Duplicated Ensembl IDs after version removal."
  )
}

background_original <- sort(unique(background_original))

cat("\n========== ORIGINAL BACKGROUND ==========\n")
cat(
  "Unique Ensembl genes tested:",
  length(background_original),
  "\n"
)


# 5. LOAD AND VALIDATE RAW MAPPING
# ============================================================

mapping_raw_path <- file.path(
  previous_mapping_dir,
  "ensembl_entrez_mapping_raw.csv"
)

if (!file.exists(mapping_raw_path)) {
  stop(
    "Previous raw mapping file not found: ",
    mapping_raw_path
  )
}

mapping_raw <- read.csv(
  mapping_raw_path,
  stringsAsFactors = FALSE
)

required_mapping_columns <- c(
  "ENSEMBL",
  "ENTREZID"
)

if (!all(required_mapping_columns %in% names(mapping_raw))) {
  stop("Invalid columns in raw mapping file.")
}

mapping_raw$ENSEMBL <- as.character(
  mapping_raw$ENSEMBL
)

mapping_raw$ENTREZID <- as.character(
  mapping_raw$ENTREZID
)

mapping_raw <- mapping_raw[
  !is.na(mapping_raw$ENSEMBL) &
    !is.na(mapping_raw$ENTREZID) &
    mapping_raw$ENSEMBL != "" &
    mapping_raw$ENTREZID != "",
]

mapping_raw <- unique(mapping_raw)

# Keep only mappings belonging to the current background.
mapping_raw <- mapping_raw[
  mapping_raw$ENSEMBL %in% background_original,
]

if (nrow(mapping_raw) == 0) {
  stop("No mapping associations found for the background.")
}


# 6. STRICT ONE-TO-ONE MAPPING
# ============================================================

ensembl_counts <- table(mapping_raw$ENSEMBL)
entrez_counts <- table(mapping_raw$ENTREZID)

unique_ensembl <- names(
  ensembl_counts[ensembl_counts == 1]
)

unique_entrez <- names(
  entrez_counts[entrez_counts == 1]
)

strict_mapping <- mapping_raw[
  mapping_raw$ENSEMBL %in% unique_ensembl &
    mapping_raw$ENTREZID %in% unique_entrez,
]

strict_mapping <- unique(strict_mapping)

strict_mapping <- strict_mapping[
  order(strict_mapping$ENSEMBL),
]

rownames(strict_mapping) <- NULL

if (anyDuplicated(strict_mapping$ENSEMBL)) {
  stop("Duplicate Ensembl IDs in strict mapping.")
}

if (anyDuplicated(strict_mapping$ENTREZID)) {
  stop("Duplicate Entrez IDs in strict mapping.")
}

if (nrow(strict_mapping) == 0) {
  stop("Strict mapping contains no usable associations.")
}

background_ensembl <- intersect(
  background_original,
  strict_mapping$ENSEMBL
)

background_mapping <- strict_mapping[
  strict_mapping$ENSEMBL %in% background_ensembl,
]

background_entrez <- unique(
  as.character(background_mapping$ENTREZID)
)

if (length(background_entrez) != nrow(background_mapping)) {
  stop("Background Entrez IDs are not unique.")
}

if (length(background_entrez) == 0) {
  stop("Strict background is empty.")
}

write.csv(
  mapping_raw,
  file.path(
    output_dir,
    "mapping",
    "ensembl_entrez_mapping_raw.csv"
  ),
  row.names = FALSE
)

write.csv(
  strict_mapping,
  file.path(
    output_dir,
    "mapping",
    "ensembl_entrez_mapping_strict.csv"
  ),
  row.names = FALSE
)

mapping_summary <- data.frame(
  metric = c(
    "Original tested genes",
    "Unique Ensembl IDs",
    "Ensembl IDs with at least one mapping",
    "Ensembl IDs with multiple Entrez mappings",
    "Entrez IDs with multiple Ensembl mappings",
    "Strict one-to-one Ensembl IDs",
    "Strict one-to-one Entrez IDs",
    "Unmapped or ambiguous Ensembl IDs",
    "Strict background coverage percentage"
  ),
  value = c(
    length(background_original),
    length(background_original),
    length(unique(mapping_raw$ENSEMBL)),
    sum(ensembl_counts > 1),
    sum(entrez_counts > 1),
    length(unique(strict_mapping$ENSEMBL)),
    length(unique(strict_mapping$ENTREZID)),
    length(setdiff(
      background_original,
      strict_mapping$ENSEMBL
    )),
    round(
      100 * length(background_ensembl) /
        length(background_original),
      2
    )
  )
)

write.csv(
  mapping_summary,
  file.path(
    output_dir,
    "mapping",
    "mapping_summary_strict.csv"
  ),
  row.names = FALSE
)

cat("\n========== STRICT MAPPING ==========\n")
print(mapping_summary, row.names = FALSE)


# 7. PREPARE SIGNIFICANT GENE LISTS
# ============================================================

prepare_gene_list <- function(data, direction, mapping) {
  
  significant <- data[
    !is.na(data$FDR) &
      data$FDR < 0.05,
  ]
  
  if (direction == "UP") {
    
    significant <- significant[
      !is.na(significant$logFC) &
        significant$logFC > 0,
    ]
    
  } else if (direction == "DOWN") {
    
    significant <- significant[
      !is.na(significant$logFC) &
        significant$logFC < 0,
    ]
    
  } else {
    stop("Unknown direction: ", direction)
  }
  
  significant$ENSEMBL <- sub(
    "\\.[0-9]+$",
    "",
    significant$gene_id
  )
  
  mapped_ensembl <- intersect(
    significant$ENSEMBL,
    mapping$ENSEMBL
  )
  
  mapped <- mapping[
    mapping$ENSEMBL %in% mapped_ensembl,
  ]
  
  entrez_ids <- unique(
    as.character(mapped$ENTREZID)
  )
  
  if (anyDuplicated(entrez_ids)) {
    stop(
      "Duplicate Entrez IDs detected in ",
      direction,
      " gene list."
    )
  }
  
  excluded_genes <- setdiff(
    significant$ENSEMBL,
    mapped_ensembl
  )
  
  list(
    genes = entrez_ids,
    n_significant = nrow(significant),
    n_mapped = length(entrez_ids),
    n_excluded = length(excluded_genes),
    coverage = if (nrow(significant) > 0) {
      100 * length(entrez_ids) / nrow(significant)
    } else {
      0
    },
    excluded_ensembl = excluded_genes
  )
}


# 8. ENRICHMENT FUNCTION
# ============================================================

run_enrichment <- function(
    gene_ids,
    background_ids,
    method
) {
  
  if (length(gene_ids) < 2) {
    return(NULL)
  }
  
  if (!all(gene_ids %in% background_ids)) {
    stop(
      "Gene list contains IDs outside the background."
    )
  }
  
  if (method %in% c("BP", "MF", "CC")) {
    
    result <- clusterProfiler::enrichGO(
      gene = gene_ids,
      universe = background_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = method,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 10,
      maxGSSize = 500,
      readable = TRUE
    )
    
  } else if (method == "REACTOME") {
    
    result <- ReactomePA::enrichPathway(
      gene = gene_ids,
      universe = background_ids,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      minGSSize = 10,
      maxGSSize = 500,
      readable = TRUE
    )
    
  } else {
    stop("Unknown enrichment method: ", method)
  }
  
  result
}


# 9. RUN ENRICHMENT ANALYSES
# ============================================================

methods <- c(
  "BP",
  "MF",
  "CC",
  "REACTOME"
)

directions <- c(
  "UP",
  "DOWN"
)

enrichment_results <- list()
summary_results <- list()
gene_coverage_results <- list()

analysis_counter <- 0

for (comparison in names(de_data)) {
  
  current_data <- de_data[[comparison]]
  
  for (direction in directions) {
    
    gene_list <- prepare_gene_list(
      data = current_data,
      direction = direction,
      mapping = strict_mapping
    )
    
    cat("\n============================================\n")
    cat("COMPARISON:", comparison, "\n")
    cat("DIRECTION:", direction, "\n")
    cat("Significant genes:", gene_list$n_significant, "\n")
    cat("Mapped genes:", gene_list$n_mapped, "\n")
    cat("Excluded genes:", gene_list$n_excluded, "\n")
    cat("Mapping coverage:", round(gene_list$coverage, 2), "%\n")
    cat("============================================\n")
    
    coverage_name <- paste(
      comparison,
      direction,
      sep = "_"
    )
    
    gene_coverage_results[[coverage_name]] <- data.frame(
      comparison = comparison,
      direction = direction,
      significant_genes = gene_list$n_significant,
      mapped_genes = gene_list$n_mapped,
      excluded_genes = gene_list$n_excluded,
      mapping_coverage_percent = round(
        gene_list$coverage,
        2
      )
    )
    
    write.csv(
      data.frame(ENSEMBL = gene_list$excluded_ensembl),
      file.path(
        output_dir,
        "reports",
        paste0(
          coverage_name,
          "_excluded_ensembl.csv"
        )
      ),
      row.names = FALSE
    )
    
    for (method in methods) {
      
      analysis_counter <- analysis_counter + 1
      
      analysis_name <- paste(
        comparison,
        direction,
        method,
        sep = "_"
      )
      
      cat("\nRunning:", analysis_name, "\n")
      
      result <- tryCatch(
        run_enrichment(
          gene_ids = gene_list$genes,
          background_ids = background_entrez,
          method = method
        ),
        error = function(e) {
          stop(
            "Enrichment failed for ",
            analysis_name,
            ": ",
            conditionMessage(e)
          )
        }
      )
      
      enrichment_results[[analysis_name]] <- result
      
      result_table <- if (is.null(result)) {
        data.frame()
      } else {
        as.data.frame(result)
      }
      
      write.csv(
        result_table,
        file.path(
          output_dir,
          "results",
          paste0(analysis_name, ".csv")
        ),
        row.names = FALSE
      )
      
      saveRDS(
        result,
        file.path(
          output_dir,
          "results",
          paste0(analysis_name, ".rds")
        )
      )
      
      n_terms <- nrow(result_table)
      
      summary_results[[analysis_name]] <- data.frame(
        comparison = comparison,
        direction = direction,
        method = method,
        significant_genes = gene_list$n_significant,
        mapped_genes = gene_list$n_mapped,
        mapping_coverage_percent = round(
          gene_list$coverage,
          2
        ),
        background_genes = length(background_entrez),
        enriched_terms = n_terms
      )
      
      cat("Enriched terms:", n_terms, "\n")
      
      if (!is.null(result) && n_terms > 0) {
        
        plot_result <- enrichplot::dotplot(
          result,
          showCategory = min(20, n_terms),
          font.size = 10
        ) +
          ggplot2::ggtitle(
            paste(
              comparison,
              direction,
              method,
              sep = " - "
            )
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "bold"
            )
          )
        
        ggsave(
          filename = file.path(
            output_dir,
            "plots",
            paste0(analysis_name, ".png")
          ),
          plot = plot_result,
          width = 11,
          height = 8,
          dpi = 300
        )
        
        ggsave(
          filename = file.path(
            output_dir,
            "plots",
            paste0(analysis_name, ".pdf")
          ),
          plot = plot_result,
          width = 11,
          height = 8
        )
      }
    }
  }
}


# 10. SAVE GLOBAL SUMMARIES
# ============================================================

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(summary_table) <- NULL

coverage_table <- do.call(
  rbind,
  gene_coverage_results
)

rownames(coverage_table) <- NULL

write.csv(
  summary_table,
  file.path(
    output_dir,
    "functional_enrichment_summary_strict.csv"
  ),
  row.names = FALSE
)

write.csv(
  coverage_table,
  file.path(
    output_dir,
    "reports",
    "gene_mapping_coverage_strict.csv"
  ),
  row.names = FALSE
)

saveRDS(
  enrichment_results,
  file.path(
    output_dir,
    "functional_enrichment_all_results_strict.rds"
  )
)


# 11. FINAL VALIDATION
# ============================================================

if (nrow(summary_table) != 16) {
  stop("Unexpected number of enrichment analyses.")
}

if (any(summary_table$background_genes !=
        length(background_entrez))) {
  stop("Inconsistent background across analyses.")
}

if (any(
  summary_table$mapped_genes >
  summary_table$significant_genes
)) {
  stop("Mapped gene count exceeds significant gene count.")
}

cat("\n============================================\n")
cat("STRICT FUNCTIONAL ENRICHMENT - FINAL SUMMARY\n")
cat("============================================\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\nGENE MAPPING COVERAGE\n")
print(
  coverage_table,
  row.names = FALSE
)

cat("\nTotal analyses:", analysis_counter, "\n")

cat(
  "Analyses with enriched terms:",
  sum(summary_table$enriched_terms > 0),
  "\n"
)

cat(
  "Analyses without enriched terms:",
  sum(summary_table$enriched_terms == 0),
  "\n"
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat("\nSTRICT FUNCTIONAL ENRICHMENT COMPLETED.\n")
