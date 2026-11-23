
# ============================================================
# TRANSCRIPTOME-X
# STEP 7 - FUNCTIONAL ENRICHMENT ANALYSIS
# ============================================================
#
# Dataset: GSE266566
# Cell line: MDA-MB-231
# Biological context: NRP1 knockdown
#
# Analyses:
#   - GO Biological Process
#   - GO Molecular Function
#   - GO Cellular Component
#   - Reactome Pathways
#
# Separate analyses:
#   - Up-regulated genes
#   - Down-regulated genes
#   - siNRP1_3 vs siNT
#   - siNRP1_5 vs siNT
#
# Statistical threshold: FDR < 0.05
# Background: all genes tested by edgeR
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

output_dir <- file.path(
  project_dir,
  "data",
  "processed",
  "functional_enrichment"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_dir, "results"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_dir, "plots"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(output_dir, "mapping"),
  recursive = TRUE,
  showWarnings = FALSE
)


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

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)


# 3. INPUT FILES
# ============================================================

input_files <- c(
  siNRP1_3_vs_siNT = "siNRP1_3_vs_siNT_annotated.csv",
  siNRP1_5_vs_siNT = "siNRP1_5_vs_siNT_annotated.csv"
)

for (f in input_files) {
  if (!file.exists(file.path(input_dir, f))) {
    stop("Input file not found: ", f)
  }
}

cat("\n============================================\n")
cat("TRANSCRIPTOME-X - FUNCTIONAL ENRICHMENT\n")
cat("============================================\n")


# 4. LOAD DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

de_data <- list()

for (comparison in names(input_files)) {
  
  file_path <- file.path(
    input_dir,
    input_files[[comparison]]
  )
  
  de_data[[comparison]] <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  if (!all(required_columns %in% names(de_data[[comparison]]))) {
    stop("Missing columns in: ", input_files[[comparison]])
  }
  
  if (anyDuplicated(de_data[[comparison]]$gene_id)) {
    stop("Duplicated gene IDs in: ", input_files[[comparison]])
  }
  
  cat(
    "\nComparison:",
    comparison,
    "\nGenes tested:",
    nrow(de_data[[comparison]]),
    "\n"
  )
}


# 5. DEFINE BACKGROUND
# ============================================================
#
# The background consists of all genes tested by edgeR.
# No filtering by FDR or logFC is applied to the background.
#
# Ensembl version suffixes are removed only for mapping.
# Original IDs are preserved in the input files.
#
# ============================================================

background_ids <- unique(
  unlist(
    lapply(
      de_data,
      function(x) x$gene_id
    )
  )
)

# Both comparisons must use the same tested genes.
if (
  !setequal(
    de_data[[1]]$gene_id,
    de_data[[2]]$gene_id
  )
) {
  stop("The two comparisons do not have the same gene universe.")
}

background_ensembl <- unique(
  sub("\\.[0-9]+$", "", background_ids)
)

background_ensembl <- background_ensembl[
  !is.na(background_ensembl) &
    background_ensembl != ""
]

cat("\n========== BACKGROUND ==========\n")
cat("Original tested genes:", length(background_ids), "\n")
cat("Unique Ensembl IDs:", length(background_ensembl), "\n")


# 6. ENSEMBL TO ENTREZ ID MAPPING
# ============================================================
#
# Mapping is performed once on the complete background.
# Ambiguous mappings are resolved deterministically:
#   1. One Entrez ID per Ensembl ID.
#   2. One Ensembl ID per Entrez ID.
#
# This prevents duplicated genes from inflating enrichment.
#
# ============================================================

mapping_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = background_ensembl,
  keytype = "ENSEMBL",
  columns = c("ENTREZID")
)

mapping_raw <- mapping_raw[
  !is.na(mapping_raw$ENTREZID) &
    !is.na(mapping_raw$ENSEMBL),
]

mapping_raw$ENSEMBL <- as.character(mapping_raw$ENSEMBL)
mapping_raw$ENTREZID <- as.character(mapping_raw$ENTREZID)

# Remove duplicated pairs.
mapping_raw <- unique(mapping_raw)

# Resolve one-to-many Ensembl mappings.
mapping_one_to_one <- mapping_raw[
  order(
    mapping_raw$ENSEMBL,
    mapping_raw$ENTREZID
  ),
]

mapping_one_to_one <- mapping_one_to_one[
  !duplicated(mapping_one_to_one$ENSEMBL),
]

# Resolve multiple Ensembl IDs mapping to the same Entrez ID.
mapping_one_to_one <- mapping_one_to_one[
  order(
    mapping_one_to_one$ENTREZID,
    mapping_one_to_one$ENSEMBL
  ),
]

mapping_one_to_one <- mapping_one_to_one[
  !duplicated(mapping_one_to_one$ENTREZID),
]

rownames(mapping_one_to_one) <- NULL

background_entrez <- unique(
  mapping_one_to_one$ENTREZID
)

if (length(background_entrez) == 0) {
  stop("No background genes could be mapped to Entrez IDs.")
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
  mapping_one_to_one,
  file.path(
    output_dir,
    "mapping",
    "ensembl_entrez_mapping_used.csv"
  ),
  row.names = FALSE
)

mapping_summary <- data.frame(
  metric = c(
    "Original background genes",
    "Unique Ensembl IDs",
    "Ensembl IDs with at least one mapping",
    "Unique Entrez IDs used",
    "Unmapped Ensembl IDs",
    "Ensembl IDs with multiple Entrez mappings"
  ),
  value = c(
    length(background_ids),
    length(background_ensembl),
    length(unique(mapping_raw$ENSEMBL)),
    length(background_entrez),
    length(setdiff(
      background_ensembl,
      unique(mapping_raw$ENSEMBL)
    )),
    sum(
      table(mapping_raw$ENSEMBL) > 1
    )
  )
)

write.csv(
  mapping_summary,
  file.path(
    output_dir,
    "mapping",
    "mapping_summary.csv"
  ),
  row.names = FALSE
)

cat("\n========== MAPPING SUMMARY ==========\n")
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
      significant$logFC > 0,
    ]
  } else if (direction == "DOWN") {
    significant <- significant[
      significant$logFC < 0,
    ]
  } else {
    stop("Unknown direction.")
  }
  
  significant$ENSEMBL <- sub(
    "\\.[0-9]+$",
    "",
    significant$gene_id
  )
  
  mapped <- merge(
    significant,
    mapping,
    by = "ENSEMBL",
    all = FALSE
  )
  
  entrez_ids <- unique(
    as.character(mapped$ENTREZID)
  )
  
  entrez_ids <- intersect(
    entrez_ids,
    mapping$ENTREZID
  )
  
  list(
    genes = entrez_ids,
    n_significant = nrow(significant),
    n_mapped = length(entrez_ids)
  )
}


# 8. ENRICHMENT FUNCTION
# ============================================================

run_enrichment <- function(
    gene_ids,
    background_ids,
    method,
    comparison,
    direction
) {
  
  if (length(gene_ids) < 2) {
    cat(
      "\nInsufficient genes for:",
      comparison,
      direction,
      method,
      "\n"
    )
    
    return(NULL)
  }
  
  if (method %in% c("BP", "MF", "CC")) {
    
    result <- tryCatch(
      clusterProfiler::enrichGO(
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
      ),
      error = function(e) {
        warning(
          "GO enrichment failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
  } else if (method == "REACTOME") {
    
    result <- tryCatch(
      ReactomePA::enrichPathway(
        gene = gene_ids,
        universe = background_ids,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      ),
      error = function(e) {
        warning(
          "Reactome enrichment failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
  } else {
    stop("Unknown enrichment method.")
  }
  
  result
}


# 9. RUN ALL ENRICHMENT ANALYSES
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

analysis_counter <- 0

for (comparison in names(de_data)) {
  
  current_data <- de_data[[comparison]]
  
  for (direction in directions) {
    
    gene_list <- prepare_gene_list(
      data = current_data,
      direction = direction,
      mapping = mapping_one_to_one
    )
    
    cat("\n============================================\n")
    cat("COMPARISON:", comparison, "\n")
    cat("DIRECTION:", direction, "\n")
    cat("Significant genes:", gene_list$n_significant, "\n")
    cat("Mapped genes:", gene_list$n_mapped, "\n")
    cat("============================================\n")
    
    for (method in methods) {
      
      analysis_counter <- analysis_counter + 1
      
      analysis_name <- paste(
        comparison,
        direction,
        method,
        sep = "_"
      )
      
      cat("\nRunning:", analysis_name, "\n")
      
      result <- run_enrichment(
        gene_ids = gene_list$genes,
        background_ids = background_entrez,
        method = method,
        comparison = comparison,
        direction = direction
      )
      
      enrichment_results[[analysis_name]] <- result
      
      result_table <- if (is.null(result)) {
        data.frame()
      } else {
        as.data.frame(result)
      }
      
      output_csv <- file.path(
        output_dir,
        "results",
        paste0(analysis_name, ".csv")
      )
      
      write.csv(
        result_table,
        output_csv,
        row.names = FALSE
      )
      
      # Save the enrichment object.
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
        background_genes = length(background_entrez),
        enriched_terms = n_terms
      )
      
      cat("Enriched terms:", n_terms, "\n")
      
      # Generate dotplot only when terms are available.
      if (!is.null(result) && n_terms > 0) {
        
        plot_result <- tryCatch(
          enrichplot::dotplot(
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
            ),
          error = function(e) {
            warning(
              "Plot generation failed for ",
              analysis_name,
              ": ",
              conditionMessage(e)
            )
            NULL
          }
        )
        
        if (!is.null(plot_result)) {
          
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
}


# 10. SAVE GLOBAL SUMMARY
# ============================================================

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(summary_table) <- NULL

write.csv(
  summary_table,
  file.path(
    output_dir,
    "functional_enrichment_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  enrichment_results,
  file.path(
    output_dir,
    "functional_enrichment_all_results.rds"
  )
)


# 11. FINAL VALIDATION
# ============================================================

cat("\n============================================\n")
cat("FUNCTIONAL ENRICHMENT - FINAL SUMMARY\n")
cat("============================================\n")

print(
  summary_table,
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

cat("\nSTEP 7 - FUNCTIONAL ENRICHMENT COMPLETED.\n")
