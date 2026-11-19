
# ============================================================
# TRANSCRIPTOME-X
# 18 - CORE GENE ANALYSIS
# ============================================================
#
# Dataset: GSE266566
# Cell model: MDA-MB-231
#
# Objectives:
# - Identify genes recurring across enriched terms
# - Analyze gene recurrence by comparison and direction
# - Separate GO and Reactome contributions
# - Integrate recurrence with differential expression
#
# Uses ALL significant enrichment terms, not only top 10.
# Does not modify previous analysis outputs.
# ============================================================


cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE ANALYSIS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

project_root <- getwd()

input_dir <- file.path(
  project_root,
  "data/processed/functional_enrichment_strict/reports",
  "integrated_biological_interpretation"
)

de_dir <- file.path(
  project_root,
  "data/processed/differential_expression/annotated"
)

output_dir <- file.path(
  project_root,
  "data/processed/functional_enrichment_strict/reports",
  "core_gene_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


comparisons <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

directions <- c("UP", "DOWN")

methods <- c("BP", "MF", "CC", "REACTOME")


# ------------------------------------------------------------
# 2. LOAD ENRICHMENT RESULTS
# ------------------------------------------------------------

cat("1. IMPORTAZIONE ENRICHMENT RESULTS\n")
cat("------------------------------------------------------------\n")


enrichment_file <- file.path(
  input_dir,
  "integrated_enrichment_results.csv"
)

if (!file.exists(enrichment_file)) {
  stop("File enrichment non trovato: ", enrichment_file)
}

enrichment <- read.csv(
  enrichment_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_enrichment_columns <- c(
  "ID",
  "Description",
  "p.adjust",
  "geneID",
  "comparison",
  "direction",
  "method"
)

missing_columns <- setdiff(
  required_enrichment_columns,
  colnames(enrichment)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti nel file enrichment: ",
    paste(missing_columns, collapse = ", ")
  )
}


enrichment <- enrichment[
  !is.na(enrichment$geneID) &
    enrichment$geneID != "",
  ,
  drop = FALSE
]


cat("Termini con geni associati:", nrow(enrichment), "\n")
cat(
  "Confronti:",
  length(unique(enrichment$comparison)),
  "\n"
)


# ------------------------------------------------------------
# 3. EXPAND GENE MEMBERSHIP
# ------------------------------------------------------------

cat("\n2. ESPANSIONE GENE MEMBERSHIP\n")
cat("------------------------------------------------------------\n")


gene_rows <- vector("list", nrow(enrichment))

for (i in seq_len(nrow(enrichment))) {
  
  genes <- unique(
    trimws(
      unlist(
        strsplit(
          enrichment$geneID[i],
          "/",
          fixed = TRUE
        )
      )
    )
  )
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  if (length(genes) == 0) {
    next
  }
  
  gene_rows[[i]] <- data.frame(
    comparison = enrichment$comparison[i],
    direction = enrichment$direction[i],
    method = enrichment$method[i],
    term_id = enrichment$ID[i],
    term_description = enrichment$Description[i],
    adjusted_p_value = enrichment$`p.adjust`[i],
    gene_name = genes,
    stringsAsFactors = FALSE
  )
}


gene_rows <- gene_rows[
  !vapply(gene_rows, is.null, logical(1))
]

if (length(gene_rows) == 0) {
  stop("Nessun gene estratto dai risultati di enrichment.")
}


gene_membership <- do.call(
  rbind,
  gene_rows
)

rownames(gene_membership) <- NULL


# Remove duplicate gene-term associations.
gene_membership <- unique(gene_membership)


write.csv(
  gene_membership,
  file.path(output_dir, "gene_pathway_membership.csv"),
  row.names = FALSE
)


cat("Associazioni gene-term:", nrow(gene_membership), "\n")
cat(
  "Geni distinti:",
  length(unique(gene_membership$gene_name)),
  "\n"
)


# ------------------------------------------------------------
# 4. GENE FREQUENCY BY COMPARISON, DIRECTION AND METHOD
# ------------------------------------------------------------

cat("\n3. CALCOLO FREQUENZA DEI CORE GENES\n")
cat("------------------------------------------------------------\n")


group_columns <- c(
  "comparison",
  "direction",
  "method"
)

group_key <- interaction(
  gene_membership$comparison,
  gene_membership$direction,
  gene_membership$method,
  drop = TRUE,
  lex.order = TRUE
)


frequency_list <- split(
  gene_membership,
  group_key
)


frequency_results <- lapply(
  frequency_list,
  function(current) {
    
    gene_groups <- split(
      current,
      current$gene_name
    )
    
    do.call(
      rbind,
      lapply(
        gene_groups,
        function(g) {
          
          data.frame(
            comparison = g$comparison[1],
            direction = g$direction[1],
            method = g$method[1],
            gene_name = g$gene_name[1],
            pathway_count = length(unique(g$term_id)),
            term_ids = paste(
              unique(g$term_id),
              collapse = ";"
            ),
            term_descriptions = paste(
              unique(g$term_description),
              collapse = " | "
            ),
            best_adjusted_p_value = min(
              g$adjusted_p_value,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }
)


gene_frequency <- do.call(
  rbind,
  frequency_results
)

rownames(gene_frequency) <- NULL


gene_frequency <- gene_frequency[
  order(
    gene_frequency$comparison,
    gene_frequency$direction,
    gene_frequency$method,
    -gene_frequency$pathway_count,
    gene_frequency$best_adjusted_p_value
  ),
]


write.csv(
  gene_frequency,
  file.path(output_dir, "core_gene_frequency_by_method.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. OVERALL FREQUENCY BY COMPARISON AND DIRECTION
# ------------------------------------------------------------

cat("\n4. FREQUENZA COMPLESSIVA PER CONFRONTO E DIREZIONE\n")
cat("------------------------------------------------------------\n")


overall_key <- interaction(
  gene_membership$comparison,
  gene_membership$direction,
  drop = TRUE,
  lex.order = TRUE
)


overall_groups <- split(
  gene_membership,
  overall_key
)


overall_results <- lapply(
  overall_groups,
  function(current) {
    
    gene_groups <- split(
      current,
      current$gene_name
    )
    
    do.call(
      rbind,
      lapply(
        gene_groups,
        function(g) {
          
          # Keep separate term IDs and method labels.
          unique_terms <- unique(
            g[, c("method", "term_id")]
          )
          
          data.frame(
            comparison = g$comparison[1],
            direction = g$direction[1],
            gene_name = g$gene_name[1],
            total_term_associations = nrow(g),
            unique_terms = nrow(unique_terms),
            methods = paste(
              sort(unique(g$method)),
              collapse = ";"
            ),
            term_ids = paste(
              paste(
                unique_terms$method,
                unique_terms$term_id,
                sep = ":"
              ),
              collapse = ";"
            ),
            best_adjusted_p_value = min(
              g$adjusted_p_value,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }
)


overall_frequency <- do.call(
  rbind,
  overall_results
)

rownames(overall_frequency) <- NULL


overall_frequency <- overall_frequency[
  order(
    overall_frequency$comparison,
    overall_frequency$direction,
    -overall_frequency$unique_terms,
    overall_frequency$best_adjusted_p_value
  ),
]


write.csv(
  overall_frequency,
  file.path(output_dir, "core_gene_frequency_overall.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. LOAD DIFFERENTIAL EXPRESSION RESULTS
# ------------------------------------------------------------

cat("\n5. IMPORTAZIONE DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")


de_list <- list()

for (comparison in comparisons) {
  
  de_file <- file.path(
    de_dir,
    paste0(comparison, "_annotated.csv")
  )
  
  if (!file.exists(de_file)) {
    stop("File DE non trovato: ", de_file)
  }
  
  de <- read.csv(
    de_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_de_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  missing_de_columns <- setdiff(
    required_de_columns,
    colnames(de)
  )
  
  if (length(missing_de_columns) > 0) {
    stop(
      "Colonne mancanti in ",
      de_file,
      ": ",
      paste(missing_de_columns, collapse = ", ")
    )
  }
  
  de$comparison <- comparison
  
  de_list[[comparison]] <- de
}


de_all <- do.call(
  rbind,
  de_list
)

de_all <- de_all[
  !is.na(de_all$gene_name) &
    de_all$gene_name != "",
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. INTEGRATE CORE GENES WITH DE RESULTS
# ------------------------------------------------------------

cat("\n6. INTEGRAZIONE CON DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")


de_mapping <- unique(
  de_all[
    ,
    c(
      "comparison",
      "gene_id",
      "gene_name",
      "logFC",
      "FDR"
    )
  ]
)


integrated_core_genes <- merge(
  overall_frequency,
  de_mapping,
  by = c("comparison", "gene_name"),
  all.x = TRUE,
  sort = FALSE
)


integrated_core_genes$de_significant <- (
  !is.na(integrated_core_genes$FDR) &
    integrated_core_genes$FDR < 0.05
)


integrated_core_genes$de_direction <- ifelse(
  is.na(integrated_core_genes$logFC),
  NA_character_,
  ifelse(
    integrated_core_genes$logFC > 0,
    "UP",
    ifelse(
      integrated_core_genes$logFC < 0,
      "DOWN",
      "ZERO"
    )
  )
)


integrated_core_genes$direction_consistent <- (
  integrated_core_genes$de_significant &
    integrated_core_genes$direction ==
    integrated_core_genes$de_direction
)


integrated_core_genes <- integrated_core_genes[
  order(
    integrated_core_genes$comparison,
    integrated_core_genes$direction,
    -integrated_core_genes$unique_terms,
    integrated_core_genes$FDR
  ),
]


write.csv(
  integrated_core_genes,
  file.path(output_dir, "core_genes_with_differential_expression.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. TOP 20 CORE GENES PER COMPARISON AND DIRECTION
# ------------------------------------------------------------

cat("\n7. IDENTIFICAZIONE TOP CORE GENES\n")
cat("------------------------------------------------------------\n")


top_core_list <- list()

for (comparison in comparisons) {
  
  for (direction in directions) {
    
    current <- integrated_core_genes[
      integrated_core_genes$comparison == comparison &
        integrated_core_genes$direction == direction,
      ,
      drop = FALSE
    ]
    
    current <- current[
      order(
        -current$unique_terms,
        current$FDR
      ),
      ,
      drop = FALSE
    ]
    
    if (nrow(current) > 0) {
      
      top_core_list[[
        paste(comparison, direction, sep = "_")
      ]] <- head(current, 20)
    }
  }
}


if (length(top_core_list) > 0) {
  
  top_core_genes <- do.call(
    rbind,
    top_core_list
  )
  
  rownames(top_core_genes) <- NULL
  
} else {
  
  top_core_genes <- integrated_core_genes[0, , drop = FALSE]
}


write.csv(
  top_core_genes,
  file.path(output_dir, "top_20_core_genes.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. MAPPING AMBIGUITY AUDIT
# ------------------------------------------------------------

cat("\n8. CONTROLLO AMBIGUITA' DEI SIMBOLI\n")
cat("------------------------------------------------------------\n")


symbol_mapping <- unique(
  de_mapping[
    ,
    c("comparison", "gene_name", "gene_id")
  ]
)


symbol_counts <- aggregate(
  gene_id ~ comparison + gene_name,
  data = symbol_mapping,
  FUN = function(x) length(unique(x))
)


names(symbol_counts)[
  names(symbol_counts) == "gene_id"
] <- "gene_id_count"


ambiguous_symbols <- symbol_counts[
  symbol_counts$gene_id_count > 1,
  ,
  drop = FALSE
]


write.csv(
  ambiguous_symbols,
  file.path(output_dir, "ambiguous_gene_symbols.csv"),
  row.names = FALSE
)


cat(
  "Associazioni simbolo-gene_id multiple:",
  nrow(ambiguous_symbols),
  "\n"
)


# ------------------------------------------------------------
# 10. SUMMARY
# ------------------------------------------------------------

cat("\n9. RIEPILOGO CORE GENE ANALYSIS\n")
cat("------------------------------------------------------------\n")


summary_table <- do.call(
  rbind,
  lapply(
    comparisons,
    function(comparison) {
      
      do.call(
        rbind,
        lapply(
          directions,
          function(direction) {
            
            current <- integrated_core_genes[
              integrated_core_genes$comparison == comparison &
                integrated_core_genes$direction == direction,
              ,
              drop = FALSE
            ]
            
            data.frame(
              comparison = comparison,
              direction = direction,
              core_genes = length(unique(current$gene_name)),
              significant_DE_genes = length(
                unique(
                  current$gene_name[
                    current$de_significant
                  ]
                )
              ),
              direction_consistent_genes = length(
                unique(
                  current$gene_name[
                    current$direction_consistent
                  ]
                )
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)


write.csv(
  summary_table,
  file.path(output_dir, "core_gene_summary.csv"),
  row.names = FALSE
)


print(summary_table, row.names = FALSE)


# ------------------------------------------------------------
# 11. MARKDOWN REPORT
# ------------------------------------------------------------

cat("\n10. GENERAZIONE REPORT\n")
cat("------------------------------------------------------------\n")


report <- c(
  "# TRANSCRIPTOME-X",
  "",
  "# Core gene analysis",
  "",
  "## 1. Methodology",
  "",
  "All significant GO and Reactome enrichment terms were used.",
  "Genes were extracted from the readable geneID field.",
  "Gene recurrence was calculated using unique enrichment term IDs.",
  "Core genes were integrated with differential expression results.",
  "",
  "## 2. Summary",
  "",
  paste(
    capture.output(
      print(summary_table, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 3. Top 20 core genes",
  ""
)


for (comparison in comparisons) {
  
  for (direction in directions) {
    
    current <- top_core_genes[
      top_core_genes$comparison == comparison &
        top_core_genes$direction == direction,
      ,
      drop = FALSE
    ]
    
    report <- c(
      report,
      paste0(
        "### ",
        comparison,
        " - ",
        direction
      ),
      ""
    )
    
    if (nrow(current) == 0) {
      
      report <- c(
        report,
        "No core genes identified.",
        ""
      )
      
    } else {
      
      report_table <- current[
        ,
        intersect(
          c(
            "gene_name",
            "unique_terms",
            "methods",
            "logFC",
            "FDR",
            "de_direction",
            "direction_consistent"
          ),
          colnames(current)
        ),
        drop = FALSE
      ]
      
      report <- c(
        report,
        paste(
          capture.output(
            print(report_table, row.names = FALSE)
          ),
          collapse = "\n"
        ),
        ""
      )
    }
  }
}


report <- c(
  report,
  "## 4. Important limitations",
  "",
  "- Gene recurrence is not an independent measure of biological importance.",
  "- GO terms are often hierarchical and redundant.",
  "- Reactome pathways may share many genes.",
  "- A gene occurring in many related terms may reflect annotation structure.",
  "- Readable gene symbols can be ambiguous; ambiguous symbol-to-Ensembl mappings are reported separately.",
  "- The same dataset is used for differential expression and enrichment.",
  "- Recurrence does not establish causality or pathway activation.",
  "- Expression direction is checked against the original differential expression results.",
  "",
  "## 5. Output files",
  "",
  "- gene_pathway_membership.csv",
  "- core_gene_frequency_by_method.csv",
  "- core_gene_frequency_overall.csv",
  "- core_genes_with_differential_expression.csv",
  "- top_20_core_genes.csv",
  "- ambiguous_gene_symbols.csv",
  "- core_gene_summary.csv"
)


report_file <- file.path(
  output_dir,
  "TRANSCRIPTOME_X_CORE_GENE_REPORT.md"
)


writeLines(
  report,
  report_file
)


# ------------------------------------------------------------
# 12. OUTPUT VALIDATION
# ------------------------------------------------------------

cat("\n11. OUTPUT VALIDATION\n")
cat("------------------------------------------------------------\n")


output_files <- c(
  "gene_pathway_membership.csv",
  "core_gene_frequency_by_method.csv",
  "core_gene_frequency_overall.csv",
  "core_genes_with_differential_expression.csv",
  "top_20_core_genes.csv",
  "ambiguous_gene_symbols.csv",
  "core_gene_summary.csv",
  "TRANSCRIPTOME_X_CORE_GENE_REPORT.md"
)


output_paths <- file.path(
  output_dir,
  output_files
)


validation <- data.frame(
  file = output_files,
  exists = file.exists(output_paths),
  size_kb = round(
    file.info(output_paths)$size / 1024,
    2
  )
)


print(validation, row.names = FALSE)


if (!all(validation$exists)) {
  stop("Validazione output non superata.")
}


if (nrow(gene_membership) == 0) {
  stop("La tabella gene-pathway membership è vuota.")
}


if (nrow(integrated_core_genes) == 0) {
  stop("Nessun core gene integrato con i risultati DE.")
}


cat("\n============================================================\n")
cat("CORE GENE ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "\nGeni distinti nei pathway:",
  length(unique(gene_membership$gene_name)),
  "\n"
)

cat(
  "Associazioni gene-term:",
  nrow(gene_membership),
  "\n"
)

cat(
  "Core genes integrati con DE:",
  nrow(integrated_core_genes),
  "\n"
)

cat("\nCartella report:\n")
cat(output_dir, "\n")

cat("\nCORE GENE ANALYSIS VALIDATED.\n")
