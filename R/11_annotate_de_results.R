
# ============================================================
# TRANSCRIPTOME-X
# STEP 4 - ANNOTATION AND VALIDATION OF DE RESULTS
# ============================================================

cat("\nTRANSCRIPTOME-X")
cat("\nANNOTATION AND VALIDATION OF DIFFERENTIAL EXPRESSION")
cat("\n====================================================\n")

# 1. PACCHETTI
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
})

# 2. PERCORSI
# ------------------------------------------------------------

base_dir <- "data/processed"

normalized_file <- file.path(
  base_dir,
  "normalized",
  "GSE266566_edgeR_TMM.rds"
)

de_dir <- file.path(
  base_dir,
  "differential_expression"
)

output_dir <- file.path(
  de_dir,
  "annotated"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. IMPORTAZIONE ANNOTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE ANNOTAZIONE\n")
cat("----------------------------\n")

if (!file.exists(normalized_file)) {
  stop("File normalizzato non trovato.")
}

dge <- readRDS(normalized_file)

annotation <- dge$genes

required_columns <- c(
  "Ensembl.114.Gene.ID",
  "Ensembl.114.Gene.Name",
  "Ensembl.114.Gene.Type"
)

if (!all(required_columns %in% colnames(annotation))) {
  stop("Le colonne di annotazione attese non sono disponibili.")
}

if (anyDuplicated(annotation$Ensembl.114.Gene.ID) > 0) {
  stop("Identificativi Ensembl duplicati nell'annotazione.")
}

cat("Geni nell'annotazione:", nrow(annotation), "\n")
cat("Colonne:", ncol(annotation), "\n")

# 4. FUNZIONE DI ANNOTAZIONE E VALIDAZIONE
# ------------------------------------------------------------

annotate_results <- function(input_file, comparison) {
  
  cat("\nConfronto:", comparison, "\n")
  
  if (!file.exists(input_file)) {
    stop("File dei risultati non trovato: ", input_file)
  }
  
  results <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_result_columns <- c(
    "gene_id",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR"
  )
  
  if (!all(required_result_columns %in% colnames(results))) {
    stop("Colonne statistiche mancanti in ", comparison)
  }
  
  if (anyDuplicated(results$gene_id) > 0) {
    stop("Gene ID duplicati nei risultati: ", comparison)
  }
  
  # Associazione esatta tra Ensembl ID e annotazione
  annotation_index <- match(
    results$gene_id,
    annotation$Ensembl.114.Gene.ID
  )
  
  unmatched <- sum(is.na(annotation_index))
  
  cat("Geni nei risultati:", nrow(results), "\n")
  cat("Geni non annotati:", unmatched, "\n")
  
  if (unmatched > 0) {
    stop(
      "Annotazione incompleta: ",
      unmatched,
      " identificativi non trovati."
    )
  }
  
  # Aggiunta annotazione senza alterare i risultati statistici
  results$gene_name <- annotation$Ensembl.114.Gene.Name[
    annotation_index
  ]
  
  results$gene_type <- annotation$Ensembl.114.Gene.Type[
    annotation_index
  ]
  
  # Riordino colonne
  results <- results[
    ,
    c(
      "gene_id",
      "gene_name",
      "gene_type",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR"
    )
  ]
  
  # Ordinamento per FDR crescente
  results <- results[
    order(results$FDR, results$PValue),
  ]
  
  # Controlli numerici
  if (any(!is.finite(results$logFC))) {
    stop("Valori logFC non finiti in ", comparison)
  }
  
  if (any(!is.finite(results$FDR))) {
    stop("Valori FDR non finiti in ", comparison)
  }
  
  if (any(results$FDR < 0 | results$FDR > 1)) {
    stop("Valori FDR fuori dall'intervallo [0,1].")
  }
  
  # Definizione dei gruppi di geni
  significant <- results[
    !is.na(results$FDR) & results$FDR < 0.05,
  ]
  
  strict <- results[
    !is.na(results$FDR) &
      results$FDR < 0.05 &
      abs(results$logFC) >= 1,
  ]
  
  upregulated <- significant[
    significant$logFC > 0,
  ]
  
  downregulated <- significant[
    significant$logFC < 0,
  ]
  
  # Salvataggio tabelle
  write.csv(
    results,
    file.path(
      output_dir,
      paste0(comparison, "_annotated.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    significant,
    file.path(
      output_dir,
      paste0(comparison, "_significant.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    strict,
    file.path(
      output_dir,
      paste0(comparison, "_strict.csv")
    ),
    row.names = FALSE
  )
  
  # Riepilogo
  summary_row <- data.frame(
    comparison = comparison,
    genes_tested = nrow(results),
    genes_annotated = sum(!is.na(results$gene_name)),
    genes_FDR_below_005 = nrow(significant),
    upregulated_FDR_below_005 = nrow(upregulated),
    downregulated_FDR_below_005 = nrow(downregulated),
    genes_FDR_below_005_logFC_at_least_1 = nrow(strict),
    stringsAsFactors = FALSE
  )
  
  cat("Geni significativi:", nrow(significant), "\n")
  cat("Up-regulated:", nrow(upregulated), "\n")
  cat("Down-regulated:", nrow(downregulated), "\n")
  cat("Geni con |logFC| >= 1:", nrow(strict), "\n")
  
  return(list(
    results = results,
    significant = significant,
    strict = strict,
    summary = summary_row
  ))
}

# 5. ANNOTAZIONE DEI DUE CONFRONTI
# ------------------------------------------------------------

cat("\n2. ANNOTAZIONE DEI RISULTATI\n")
cat("----------------------------\n")

analysis_3 <- annotate_results(
  file.path(de_dir, "siNRP1_3_vs_siNT_all_genes.csv"),
  "siNRP1_3_vs_siNT"
)

analysis_5 <- annotate_results(
  file.path(de_dir, "siNRP1_5_vs_siNT_all_genes.csv"),
  "siNRP1_5_vs_siNT"
)

# 6. RIEPILOGO GENERALE
# ------------------------------------------------------------

summary_results <- rbind(
  analysis_3$summary,
  analysis_5$summary
)

write.csv(
  summary_results,
  file.path(
    output_dir,
    "annotation_validation_summary.csv"
  ),
  row.names = FALSE
)

cat("\n3. RIEPILOGO FINALE\n")
cat("----------------------------\n")

print(summary_results)

# 7. CONTROLLO SPECIFICO NRP1
# ------------------------------------------------------------

cat("\n4. CONTROLLO SPECIFICO NRP1\n")
cat("----------------------------\n")

nrp1_results <- rbind(
  analysis_3$results[
    analysis_3$results$gene_name == "NRP1",
  ],
  analysis_5$results[
    analysis_5$results$gene_name == "NRP1",
  ]
)

if (nrow(nrp1_results) == 0) {
  stop("NRP1 non trovato nelle tabelle annotate.")
}

print(
  nrp1_results[
    ,
    c(
      "gene_id",
      "gene_name",
      "gene_type",
      "logFC",
      "PValue",
      "FDR"
    )
  ]
)

# 8. VERIFICA OUTPUT
# ------------------------------------------------------------

expected_files <- c(
  "siNRP1_3_vs_siNT_annotated.csv",
  "siNRP1_3_vs_siNT_significant.csv",
  "siNRP1_3_vs_siNT_strict.csv",
  "siNRP1_5_vs_siNT_annotated.csv",
  "siNRP1_5_vs_siNT_significant.csv",
  "siNRP1_5_vs_siNT_strict.csv",
  "annotation_validation_summary.csv"
)

file_status <- file.exists(
  file.path(output_dir, expected_files)
)

cat("\n5. VALIDAZIONE FILE\n")
cat("----------------------------\n")

print(
  data.frame(
    file = expected_files,
    exists = file_status
  )
)

if (!all(file_status)) {
  stop("Uno o più file di output non sono stati generati.")
}

cat("\n====================================================\n")
cat("STEP 4 COMPLETATO CORRETTAMENTE\n")
cat("Tabelle annotate e validate.\n")
cat("Output:", output_dir, "\n")
cat("====================================================\n")
