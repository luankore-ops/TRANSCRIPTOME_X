
# ============================================================
# TRANSCRIPTOME-X
# 17 - INTEGRATED BIOLOGICAL INTERPRETATION
# ============================================================
#
# Dataset: GSE266566
# Cell model: MDA-MB-231
#
# Comparisons:
# 1. siNRP1_3 vs siNT
# 2. siNRP1_5 vs siNT
#
# Differential expression: edgeR quasi-likelihood
# Functional enrichment: strict 1:1 Ensembl-Entrez mapping
#
# This script does not modify previous analysis outputs.
# ============================================================


cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - INTEGRATED BIOLOGICAL INTERPRETATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. PROJECT CONFIGURATION
# ------------------------------------------------------------

project_root <- getwd()

required_directories <- c(
  "data/processed/differential_expression/annotated",
  "data/processed/functional_enrichment_strict/results",
  "data/processed/functional_enrichment_strict/reports"
)

missing_directories <- required_directories[
  !dir.exists(file.path(project_root, required_directories))
]

if (length(missing_directories) > 0) {
  stop(
    "Cartelle mancanti:\n",
    paste(missing_directories, collapse = "\n")
  )
}


output_dir <- file.path(
  project_root,
  "data/processed/functional_enrichment_strict/reports",
  "integrated_biological_interpretation"
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

methods <- c(
  "BP",
  "MF",
  "CC",
  "REACTOME"
)


# ------------------------------------------------------------
# 2. LOAD DIFFERENTIAL EXPRESSION RESULTS
# ------------------------------------------------------------

cat("1. IMPORTAZIONE DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")


de_list <- list()

for (comparison in comparisons) {
  
  input_file <- file.path(
    project_root,
    "data/processed/differential_expression/annotated",
    paste0(comparison, "_annotated.csv")
  )
  
  if (!file.exists(input_file)) {
    stop("File DE non trovato: ", input_file)
  }
  
  de <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(de)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Colonne mancanti nel file ",
      comparison,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  de$comparison <- rep(comparison, nrow(de))
  
  de_list[[comparison]] <- de
  
  cat(
    comparison,
    "| Geni:",
    nrow(de),
    "| FDR < 0.05:",
    sum(de$FDR < 0.05, na.rm = TRUE),
    "\n"
  )
}


de_all <- do.call(rbind, de_list)

if (is.null(de_all) || nrow(de_all) == 0) {
  stop("Nessun risultato di differential expression disponibile.")
}


# ------------------------------------------------------------
# 3. LOAD STRICT FUNCTIONAL ENRICHMENT RESULTS
# ------------------------------------------------------------

cat("\n2. IMPORTAZIONE FUNCTIONAL ENRICHMENT\n")
cat("------------------------------------------------------------\n")


enrichment_list <- list()

for (comparison in comparisons) {
  
  for (direction in directions) {
    
    for (method in methods) {
      
      file_name <- paste0(
        comparison,
        "_",
        direction,
        "_",
        method,
        ".csv"
      )
      
      input_file <- file.path(
        project_root,
        "data/processed/functional_enrichment_strict/results",
        file_name
      )
      
      if (!file.exists(input_file)) {
        stop("Risultato enrichment non trovato: ", input_file)
      }
      
      result <- read.csv(
        input_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      expected_columns <- c(
        "ID",
        "Description",
        "p.adjust",
        "geneID"
      )
      
      missing_columns <- setdiff(
        expected_columns,
        colnames(result)
      )
      
      if (length(missing_columns) > 0) {
        stop(
          "Colonne mancanti in ",
          file_name,
          ": ",
          paste(missing_columns, collapse = ", ")
        )
      }
      
      # Gestione corretta anche delle tabelle vuote.
      result$comparison <- rep(comparison, nrow(result))
      result$direction <- rep(direction, nrow(result))
      result$method <- rep(method, nrow(result))
      
      enrichment_list[[file_name]] <- result
      
      cat(
        comparison,
        "|",
        direction,
        "|",
        method,
        "| Termini:",
        nrow(result),
        "\n"
      )
    }
  }
}


# Unisce anche le tabelle vuote senza interrompere l'esecuzione.
enrichment_all <- do.call(
  rbind,
  enrichment_list
)

if (is.null(enrichment_all)) {
  stop("Impossibile integrare i risultati di enrichment.")
}

rownames(enrichment_all) <- NULL


# ------------------------------------------------------------
# 4. CREATE INTEGRATED ENRICHMENT TABLE
# ------------------------------------------------------------

cat("\n3. INTEGRAZIONE DEI RISULTATI\n")
cat("------------------------------------------------------------\n")


enrichment_all <- enrichment_all[
  order(
    enrichment_all$comparison,
    enrichment_all$direction,
    enrichment_all$method,
    enrichment_all$`p.adjust`
  ),
]


write.csv(
  enrichment_all,
  file.path(output_dir, "integrated_enrichment_results.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. TOP ENRICHED TERMS
# ------------------------------------------------------------

cat("\n4. IDENTIFICAZIONE TOP ENRICHED TERMS\n")
cat("------------------------------------------------------------\n")


top_terms <- list()

for (comparison in comparisons) {
  
  for (direction in directions) {
    
    for (method in methods) {
      
      subset_result <- enrichment_all[
        enrichment_all$comparison == comparison &
          enrichment_all$direction == direction &
          enrichment_all$method == method,
        ,
        drop = FALSE
      ]
      
      if (nrow(subset_result) > 0) {
        
        subset_result <- subset_result[
          order(subset_result$`p.adjust`),
          ,
          drop = FALSE
        ]
        
        top_terms[[paste(
          comparison,
          direction,
          method,
          sep = "_"
        )]] <- head(subset_result, 10)
      }
    }
  }
}


if (length(top_terms) > 0) {
  
  top_terms_table <- do.call(
    rbind,
    top_terms
  )
  
  rownames(top_terms_table) <- NULL
  
} else {
  
  top_terms_table <- enrichment_all[0, , drop = FALSE]
}


write.csv(
  top_terms_table,
  file.path(output_dir, "top_enriched_terms.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. COMPARE ENRICHED TERMS BETWEEN CONTRASTS
# ------------------------------------------------------------

cat("\n5. CONFRONTO TRA I DUE KNOCKDOWN\n")
cat("------------------------------------------------------------\n")


overlap_list <- list()

for (direction in directions) {
  
  for (method in methods) {
    
    result_3 <- enrichment_all[
      enrichment_all$comparison == "siNRP1_3_vs_siNT" &
        enrichment_all$direction == direction &
        enrichment_all$method == method,
      ,
      drop = FALSE
    ]
    
    result_5 <- enrichment_all[
      enrichment_all$comparison == "siNRP1_5_vs_siNT" &
        enrichment_all$direction == direction &
        enrichment_all$method == method,
      ,
      drop = FALSE
    ]
    
    ids_3 <- unique(result_3$ID)
    ids_5 <- unique(result_5$ID)
    
    shared <- intersect(ids_3, ids_5)
    only_3 <- setdiff(ids_3, ids_5)
    only_5 <- setdiff(ids_5, ids_3)
    
    union_ids <- union(ids_3, ids_5)
    
    jaccard <- if (length(union_ids) > 0) {
      length(shared) / length(union_ids)
    } else {
      NA_real_
    }
    
    overlap_list[[paste(direction, method, sep = "_")]] <-
      data.frame(
        direction = direction,
        method = method,
        terms_siNRP1_3 = length(ids_3),
        terms_siNRP1_5 = length(ids_5),
        shared_terms = length(shared),
        only_siNRP1_3 = length(only_3),
        only_siNRP1_5 = length(only_5),
        jaccard = jaccard,
        stringsAsFactors = FALSE
      )
  }
}


overlap_table <- do.call(
  rbind,
  overlap_list
)

rownames(overlap_table) <- NULL


write.csv(
  overlap_table,
  file.path(output_dir, "cross_comparison_enrichment_overlap.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. NRP1 EXPRESSION SUMMARY
# ------------------------------------------------------------

cat("\n6. VERIFICA DEL GENE NRP1\n")
cat("------------------------------------------------------------\n")


nrp1_results <- de_all[
  de_all$gene_name == "NRP1" |
    grepl("^ENSG00000099250", de_all$gene_id),
  ,
  drop = FALSE
]


if (nrow(nrp1_results) == 0) {
  
  warning("NRP1 non trovato nei risultati annotati.")
  
  nrp1_results <- data.frame()
  
} else {
  
  nrp1_results <- nrp1_results[
    ,
    c(
      "comparison",
      "gene_id",
      "gene_name",
      "logFC",
      "FDR"
    ),
    drop = FALSE
  ]
  
  nrp1_results$fold_change <- 2 ^ nrp1_results$logFC
  
  write.csv(
    nrp1_results,
    file.path(output_dir, "NRP1_expression_summary.csv"),
    row.names = FALSE
  )
  
  print(nrp1_results)
}


if (nrow(nrp1_results) == 0) {
  
  write.csv(
    data.frame(
      comparison = character(),
      gene_id = character(),
      gene_name = character(),
      logFC = numeric(),
      FDR = numeric(),
      fold_change = numeric()
    ),
    file.path(output_dir, "NRP1_expression_summary.csv"),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 8. DIFFERENTIAL EXPRESSION SUMMARY
# ------------------------------------------------------------

cat("\n7. RIEPILOGO DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")


de_summary <- do.call(
  rbind,
  lapply(comparisons, function(comparison) {
    
    current <- de_all[
      de_all$comparison == comparison,
      ,
      drop = FALSE
    ]
    
    significant <- !is.na(current$FDR) &
      current$FDR < 0.05
    
    data.frame(
      comparison = comparison,
      tested_genes = nrow(current),
      significant_genes = sum(significant),
      up_regulated = sum(
        significant & current$logFC > 0,
        na.rm = TRUE
      ),
      down_regulated = sum(
        significant & current$logFC < 0,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
)


write.csv(
  de_summary,
  file.path(output_dir, "differential_expression_summary.csv"),
  row.names = FALSE
)


print(de_summary, row.names = FALSE)


# ------------------------------------------------------------
# 9. CREATE MARKDOWN REPORT
# ------------------------------------------------------------

cat("\n8. GENERAZIONE REPORT INTEGRATO\n")
cat("------------------------------------------------------------\n")


report_file <- file.path(
  output_dir,
  "TRANSCRIPTOME_X_INTEGRATED_BIOLOGICAL_INTERPRETATION.md"
)


report <- c(
  "# TRANSCRIPTOME-X",
  "",
  "# Integrated biological interpretation",
  "",
  "## 1. Study overview",
  "",
  "Dataset: GSE266566.",
  "",
  "Cell model: MDA-MB-231.",
  "",
  "Experimental conditions:",
  "- siNT: non-targeting control",
  "- siNRP1_3: NRP1 knockdown condition 3",
  "- siNRP1_5: NRP1 knockdown condition 5",
  "",
  "Comparisons:",
  "- siNRP1_3 vs siNT",
  "- siNRP1_5 vs siNT",
  "",
  "Differential expression was assessed using edgeR quasi-likelihood testing.",
  "Functional enrichment uses strict one-to-one Ensembl-Entrez mapping.",
  "",
  "## 2. Differential expression summary",
  "",
  paste(
    capture.output(
      print(de_summary, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 3. NRP1 expression",
  ""
)


if (nrow(nrp1_results) > 0) {
  
  report <- c(
    report,
    paste(
      capture.output(
        print(nrp1_results, row.names = FALSE)
      ),
      collapse = "\n"
    )
  )
  
} else {
  
  report <- c(
    report,
    "NRP1 results were not available."
  )
}


report <- c(
  report,
  "",
  "## 4. Cross-comparison enrichment overlap",
  "",
  paste(
    capture.output(
      print(overlap_table, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 5. Most significant enriched terms",
  ""
)


for (comparison in comparisons) {
  
  for (direction in directions) {
    
    for (method in methods) {
      
      current <- top_terms_table[
        top_terms_table$comparison == comparison &
          top_terms_table$direction == direction &
          top_terms_table$method == method,
        ,
        drop = FALSE
      ]
      
      report <- c(
        report,
        paste0(
          "### ",
          comparison,
          " - ",
          direction,
          " - ",
          method
        ),
        ""
      )
      
      if (nrow(current) == 0) {
        
        report <- c(
          report,
          "No significantly enriched terms were identified.",
          ""
        )
        
      } else {
        
        available_columns <- intersect(
          c("ID", "Description", "p.adjust", "Count"),
          colnames(current)
        )
        
        current_table <- current[
          ,
          available_columns,
          drop = FALSE
        ]
        
        report <- c(
          report,
          paste(
            capture.output(
              print(current_table, row.names = FALSE)
            ),
            collapse = "\n"
          ),
          ""
        )
      }
    }
  }
}


report <- c(
  report,
  "## 6. Interpretation guidelines",
  "",
  "- Enrichment indicates overrepresentation, not pathway activation or inhibition.",
  "- Differential expression does not establish a causal mechanism.",
  "- Results describe transcriptional associations in the studied cell model.",
  "- GO terms may be redundant or biologically overlapping.",
  "- The two knockdown conditions are compared separately against siNT.",
  "- Their results do not constitute a direct statistical comparison between siNRP1_3 and siNRP1_5.",
  "- The two comparisons use the same dataset and are not independent biological replications.",
  "- Enrichment results depend on annotation coverage and the selected gene background.",
  "",
  "## 7. Limitations",
  "",
  "- RNA-seq data are RSEM expected counts, including fractional values.",
  "- No RSEM inferential replicates were available to model read-assignment uncertainty.",
  "- Strict one-to-one mapping excludes unmapped and ambiguous gene identifiers.",
  "- Functional enrichment is based on genes passing the selected significance threshold.",
  "",
  "## 8. Output files",
  "",
  "- integrated_enrichment_results.csv",
  "- top_enriched_terms.csv",
  "- cross_comparison_enrichment_overlap.csv",
  "- NRP1_expression_summary.csv",
  "- differential_expression_summary.csv",
  "",
  "This report is a computational summary. Biological conclusions require critical inspection of the enriched terms, their contributing genes, and the relevant literature."
)


writeLines(
  report,
  report_file
)


# ------------------------------------------------------------
# 10. OUTPUT VALIDATION
# ------------------------------------------------------------

cat("\n9. OUTPUT VALIDATION\n")
cat("------------------------------------------------------------\n")


output_files <- c(
  "integrated_enrichment_results.csv",
  "top_enriched_terms.csv",
  "cross_comparison_enrichment_overlap.csv",
  "NRP1_expression_summary.csv",
  "differential_expression_summary.csv",
  "TRANSCRIPTOME_X_INTEGRATED_BIOLOGICAL_INTERPRETATION.md"
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
  stop("Validazione degli output non superata.")
}


if (nrow(enrichment_all) == 0) {
  stop("La tabella integrata di enrichment è vuota.")
}


# Verifica presenza di entrambi i confronti.
if (!all(comparisons %in% unique(de_all$comparison))) {
  stop("Uno o più confronti DE non sono stati importati.")
}


# Verifica che siano state importate tutte le 16 analisi.
if (length(enrichment_list) != 16) {
  stop("Numero di analisi enrichment diverso da 16.")
}


cat("\n============================================================\n")
cat("INTEGRATED BIOLOGICAL INTERPRETATION COMPLETED\n")
cat("============================================================\n")

cat(
  "\nAnalisi enrichment integrate:",
  length(enrichment_list),
  "\n"
)

cat(
  "Termini significativi complessivi:",
  nrow(enrichment_all),
  "\n"
)

cat(
  "Analisi senza termini arricchiti:",
  sum(
    vapply(
      enrichment_list,
      nrow,
      integer(1)
    ) == 0
  ),
  "\n"
)

cat(
  "Pathway confrontabili tra i due knockdown:",
  sum(!is.na(overlap_table$jaccard)),
  "\n"
)

cat("\nCartella report:\n")
cat(output_dir, "\n")

cat("\nINTEGRATED REPORT VALIDATED.\n")
