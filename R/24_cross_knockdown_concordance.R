
# ============================================================
# TRANSCRIPTOME-X - CROSS-KNOCKDOWN CONCORDANCE ANALYSIS
# Script: 24_cross_knockdown_concordance.R
#
# Obiettivi:
# 1. Confrontare tutti i geni nei due knockdown NRP1
# 2. Identificare geni significativi in entrambi i confronti
# 3. Classificare concordanza e discordanza del logFC
# 4. Identificare geni significativi in un solo confronto
# 5. Calcolare correlazioni tra i logFC
# 6. Integrare i 16 core genes concordanti identificati in R23
# 7. Generare classifiche e report
#
# Input:
# - DE annotated siNRP1_3_vs_siNT
# - DE annotated siNRP1_5_vs_siNT
# - Core genes concordanti R23
#
# Non modifica gli output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CROSS-KNOCKDOWN CONCORDANCE\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

de_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "annotated"
)

core_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports",
  "core_gene_prioritization"
)

output_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "cross_knockdown_concordance"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

file_3 <- file.path(
  de_dir,
  "siNRP1_3_vs_siNT_annotated.csv"
)

file_5 <- file.path(
  de_dir,
  "siNRP1_5_vs_siNT_annotated.csv"
)

core_file <- file.path(
  core_dir,
  "concordant_core_genes.csv"
)

required_files <- c(file_3, file_5, core_file)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "File mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Soglia statistica
FDR_THRESHOLD <- 0.05

# ------------------------------------------------------------
# 2. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")

de_3 <- read.csv(
  file_3,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

de_5 <- read.csv(
  file_5,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_16 <- read.csv(
  core_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("siNRP1_3:", nrow(de_3), "righe\n")
cat("siNRP1_5:", nrow(de_5), "righe\n")
cat("Core genes concordanti R23:", nrow(core_16), "\n")

# ------------------------------------------------------------
# 3. VALIDAZIONE COLONNE
# ------------------------------------------------------------

cat("\n2. VALIDAZIONE COLONNE\n")
cat("------------------------------------------------------------\n")

required_de_columns <- c(
  "gene_id",
  "gene_name",
  "logFC",
  "FDR"
)

for (dataset_name in c("de_3", "de_5")) {
  
  dataset <- get(dataset_name)
  
  missing <- setdiff(
    required_de_columns,
    names(dataset)
  )
  
  if (length(missing) > 0) {
    stop(
      "Colonne mancanti in ",
      dataset_name,
      ": ",
      paste(missing, collapse = ", ")
    )
  }
}

if (!"ensembl_id" %in% names(core_16)) {
  stop(
    "Colonna ensembl_id mancante nel file concordant_core_genes.csv."
  )
}

# ------------------------------------------------------------
# 4. NORMALIZZAZIONE IDENTIFICATIVI
# ------------------------------------------------------------

cat("\n3. NORMALIZZAZIONE IDENTIFICATIVI ENSEMBL\n")
cat("------------------------------------------------------------\n")

prepare_de <- function(x, comparison_name) {
  
  x$ensembl_id <- trimws(
    as.character(x$gene_id)
  )
  
  x$ensembl_id <- sub(
    "\\.[0-9]+$",
    "",
    x$ensembl_id
  )
  
  x$gene_name <- trimws(
    as.character(x$gene_name)
  )
  
  x$logFC <- as.numeric(x$logFC)
  x$FDR <- as.numeric(x$FDR)
  
  if (anyNA(x$ensembl_id) ||
      any(!nzchar(x$ensembl_id))) {
    stop(
      "Identificativi Ensembl mancanti in ",
      comparison_name
    )
  }
  
  if (anyDuplicated(x$ensembl_id) > 0) {
    stop(
      "Identificativi Ensembl duplicati in ",
      comparison_name
    )
  }
  
  if (any(!is.finite(x$logFC))) {
    stop(
      "Valori logFC non validi in ",
      comparison_name
    )
  }
  
  if (any(
    !is.na(x$FDR) &
    (x$FDR < 0 | x$FDR > 1)
  )) {
    stop(
      "Valori FDR fuori intervallo in ",
      comparison_name
    )
  }
  
  x
}

de_3 <- prepare_de(
  de_3,
  "siNRP1_3_vs_siNT"
)

de_5 <- prepare_de(
  de_5,
  "siNRP1_5_vs_siNT"
)

cat(
  "Ensembl distinti siNRP1_3:",
  length(unique(de_3$ensembl_id)),
  "\n"
)

cat(
  "Ensembl distinti siNRP1_5:",
  length(unique(de_5$ensembl_id)),
  "\n"
)

# ------------------------------------------------------------
# 5. MERGE DEI DUE CONFRONTI
# ------------------------------------------------------------

cat("\n4. INTEGRAZIONE DEI DUE KNOCKDOWN\n")
cat("------------------------------------------------------------\n")

de3_selected <- de_3[
  ,
  c(
    "ensembl_id",
    "gene_name",
    "logFC",
    "FDR"
  )
]

names(de3_selected) <- c(
  "ensembl_id",
  "gene_name_3",
  "logFC_3",
  "FDR_3"
)

de5_selected <- de_5[
  ,
  c(
    "ensembl_id",
    "gene_name",
    "logFC",
    "FDR"
  )
]

names(de5_selected) <- c(
  "ensembl_id",
  "gene_name_5",
  "logFC_5",
  "FDR_5"
)

comparison <- merge(
  de3_selected,
  de5_selected,
  by = "ensembl_id",
  all = TRUE,
  sort = FALSE
)

comparison$gene_name <- ifelse(
  !is.na(comparison$gene_name_3) &
    nzchar(comparison$gene_name_3),
  comparison$gene_name_3,
  comparison$gene_name_5
)

cat("Geni complessivi:", nrow(comparison), "\n")
cat(
  "Presenti in entrambi i confronti:",
  sum(
    !is.na(comparison$logFC_3) &
      !is.na(comparison$logFC_5)
  ),
  "\n"
)

# ------------------------------------------------------------
# 6. SIGNIFICATIVITA DIFFERENTIAL EXPRESSION
# ------------------------------------------------------------

cat("\n5. CLASSIFICAZIONE SIGNIFICATIVITA\n")
cat("------------------------------------------------------------\n")

comparison$significant_3 <- (
  !is.na(comparison$FDR_3) &
    comparison$FDR_3 < FDR_THRESHOLD
)

comparison$significant_5 <- (
  !is.na(comparison$FDR_5) &
    comparison$FDR_5 < FDR_THRESHOLD
)

comparison$direction_3 <- ifelse(
  comparison$logFC_3 > 0,
  "UP",
  ifelse(comparison$logFC_3 < 0, "DOWN", "ZERO")
)

comparison$direction_5 <- ifelse(
  comparison$logFC_5 > 0,
  "UP",
  ifelse(comparison$logFC_5 < 0, "DOWN", "ZERO")
)

comparison$significance_category <- ifelse(
  comparison$significant_3 &
    comparison$significant_5,
  "SIGNIFICANT_BOTH",
  ifelse(
    comparison$significant_3,
    "SIGNIFICANT_SI_NRP1_3_ONLY",
    ifelse(
      comparison$significant_5,
      "SIGNIFICANT_SI_NRP1_5_ONLY",
      "NOT_SIGNIFICANT"
    )
  )
)

# ------------------------------------------------------------
# 7. CONCORDANZA DELLA DIREZIONE
# ------------------------------------------------------------

cat("\n6. CONCORDANZA DELLA DIREZIONE\n")
cat("------------------------------------------------------------\n")

comparison$direction_category <- ifelse(
  comparison$direction_3 == comparison$direction_5,
  "CONCORDANT",
  "DISCORDANT"
)

comparison$cross_knockdown_status <- ifelse(
  comparison$significance_category != "SIGNIFICANT_BOTH",
  comparison$significance_category,
  ifelse(
    comparison$direction_category == "CONCORDANT",
    "BOTH_SIGNIFICANT_CONCORDANT",
    "BOTH_SIGNIFICANT_DISCORDANT"
  )
)

# La categoria concordante/discordante è biologicamente
# interpretabile come risposta replicata solo quando entrambi
# i confronti sono significativi.

comparison$replicated_direction <- (
  comparison$significance_category == "SIGNIFICANT_BOTH" &
    comparison$direction_category == "CONCORDANT" &
    comparison$direction_3 != "ZERO"
)

comparison$replicated_opposite_direction <- (
  comparison$significance_category == "SIGNIFICANT_BOTH" &
    comparison$direction_category == "DISCORDANT"
)

# ------------------------------------------------------------
# 8. RIEPILOGO CATEGORIE
# ------------------------------------------------------------

cat("\n7. DISTRIBUZIONE CATEGORIE\n")
cat("------------------------------------------------------------\n")

category_summary <- as.data.frame(
  table(comparison$cross_knockdown_status)
)

names(category_summary) <- c(
  "category",
  "gene_count"
)

category_summary$percentage <- round(
  100 * category_summary$gene_count /
    nrow(comparison),
  2
)

print(category_summary, row.names = FALSE)

write.csv(
  category_summary,
  file.path(output_dir, "cross_knockdown_category_summary.csv"),
  row.names = FALSE
)

# Tabella 2x2 di significatività
significance_table <- table(
  siNRP1_3 = comparison$significant_3,
  siNRP1_5 = comparison$significant_5
)

write.csv(
  as.data.frame.matrix(significance_table),
  file.path(output_dir, "significance_contingency_table.csv")
)

# ------------------------------------------------------------
# 9. CORRELAZIONE LOGFC
# ------------------------------------------------------------

cat("\n8. CORRELAZIONE LOGFC\n")
cat("------------------------------------------------------------\n")

valid_logfc <- (
  is.finite(comparison$logFC_3) &
    is.finite(comparison$logFC_5)
)

both_significant <- (
  valid_logfc &
    comparison$significant_3 &
    comparison$significant_5
)

concordant_significant <- (
  both_significant &
    comparison$replicated_direction
)

calculate_correlation <- function(x, y, label) {
  
  valid <- is.finite(x) & is.finite(y)
  
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 3) {
    return(
      data.frame(
        subset = label,
        n = length(x),
        pearson_r = NA_real_,
        pearson_p = NA_real_,
        spearman_rho = NA_real_,
        spearman_p = NA_real_
      )
    )
  }
  
  pearson <- cor.test(
    x,
    y,
    method = "pearson"
  )
  
  spearman <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    subset = label,
    n = length(x),
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value,
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value
  )
}

correlation_results <- rbind(
  calculate_correlation(
    comparison$logFC_3[valid_logfc],
    comparison$logFC_5[valid_logfc],
    "ALL_GENES"
  ),
  calculate_correlation(
    comparison$logFC_3[both_significant],
    comparison$logFC_5[both_significant],
    "SIGNIFICANT_BOTH"
  ),
  calculate_correlation(
    comparison$logFC_3[concordant_significant],
    comparison$logFC_5[concordant_significant],
    "SIGNIFICANT_BOTH_CONCORDANT"
  )
)

print(correlation_results, row.names = FALSE)

write.csv(
  correlation_results,
  file.path(output_dir, "logFC_correlation_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. ESTRAZIONE SOTTOINSIEMI
# ------------------------------------------------------------

cat("\n9. ESTRAZIONE SOTTOINSIEMI\n")
cat("------------------------------------------------------------\n")

significant_both <- comparison[
  comparison$significance_category == "SIGNIFICANT_BOTH",
  ,
  drop = FALSE
]

concordant <- comparison[
  comparison$replicated_direction,
  ,
  drop = FALSE
]

discordant <- comparison[
  comparison$replicated_opposite_direction,
  ,
  drop = FALSE
]

significant_3_only <- comparison[
  comparison$significance_category ==
    "SIGNIFICANT_SI_NRP1_3_ONLY",
  ,
  drop = FALSE
]

significant_5_only <- comparison[
  comparison$significance_category ==
    "SIGNIFICANT_SI_NRP1_5_ONLY",
  ,
  drop = FALSE
]

not_significant <- comparison[
  comparison$significance_category == "NOT_SIGNIFICANT",
  ,
  drop = FALSE
]

write.csv(
  comparison,
  file.path(output_dir, "all_genes_cross_knockdown.csv"),
  row.names = FALSE
)

write.csv(
  significant_both,
  file.path(output_dir, "genes_significant_both.csv"),
  row.names = FALSE
)

write.csv(
  concordant,
  file.path(output_dir, "genes_concordant_both_knockdowns.csv"),
  row.names = FALSE
)

write.csv(
  discordant,
  file.path(output_dir, "genes_discordant_both_knockdowns.csv"),
  row.names = FALSE
)

write.csv(
  significant_3_only,
  file.path(output_dir, "genes_significant_siNRP1_3_only.csv"),
  row.names = FALSE
)

write.csv(
  significant_5_only,
  file.path(output_dir, "genes_significant_siNRP1_5_only.csv"),
  row.names = FALSE
)

write.csv(
  not_significant,
  file.path(output_dir, "genes_not_significant_both.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. PRIORITIZZAZIONE GENI CONCORDANTI
# ------------------------------------------------------------

cat("\n10. PRIORITIZZAZIONE GENI CONCORDANTI\n")
cat("------------------------------------------------------------\n")

if (nrow(concordant) > 0) {
  
  concordant$minimum_FDR <- pmin(
    concordant$FDR_3,
    concordant$FDR_5
  )
  
  concordant$maximum_FDR <- pmax(
    concordant$FDR_3,
    concordant$FDR_5
  )
  
  concordant$mean_abs_logFC <- (
    abs(concordant$logFC_3) +
      abs(concordant$logFC_5)
  ) / 2
  
  concordant$minimum_abs_logFC <- pmin(
    abs(concordant$logFC_3),
    abs(concordant$logFC_5)
  )
  
  concordant$logFC_difference <- abs(
    concordant$logFC_3 -
      concordant$logFC_5
  )
  
  concordant <- concordant[
    order(
      concordant$maximum_FDR,
      -concordant$minimum_abs_logFC
    ),
  ]
  
  concordant$replication_rank <- seq_len(
    nrow(concordant)
  )
}

write.csv(
  concordant,
  file.path(output_dir, "concordant_gene_priority_ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. INTEGRAZIONE DEI 16 CORE GENES R23
# ------------------------------------------------------------

cat("\n11. INTEGRAZIONE CORE GENES R23\n")
cat("------------------------------------------------------------\n")

core_16$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  trimws(as.character(core_16$ensembl_id))
)

if (anyDuplicated(core_16$ensembl_id) > 0) {
  stop("Duplicati Ensembl nel file core concordant R23.")
}

core_cross <- merge(
  core_16,
  comparison,
  by = "ensembl_id",
  all.x = TRUE,
  sort = FALSE
)

core_cross$cross_knockdown_confirmed <- (
  !is.na(core_cross$replicated_direction) &
    core_cross$replicated_direction
)

cat(
  "Core genes R23:",
  nrow(core_16),
  "\n"
)

cat(
  "Core genes con conferma DE indipendente:",
  sum(core_cross$cross_knockdown_confirmed, na.rm = TRUE),
  "\n"
)

write.csv(
  core_cross,
  file.path(output_dir, "R23_core_genes_cross_knockdown_validation.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. GRAFICO CORRELAZIONE
# ------------------------------------------------------------

cat("\n12. GENERAZIONE GRAFICO\n")
cat("------------------------------------------------------------\n")

pdf(
  file.path(output_dir, "logFC_cross_knockdown_scatter.pdf"),
  width = 8,
  height = 7
)

plot(
  comparison$logFC_3[valid_logfc],
  comparison$logFC_5[valid_logfc],
  pch = 16,
  cex = 0.45,
  xlab = "logFC siNRP1_3 vs siNT",
  ylab = "logFC siNRP1_5 vs siNT",
  main = "Cross-knockdown logFC concordance"
)

abline(
  h = 0,
  v = 0,
  lty = 2
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

if (sum(valid_logfc) >= 3) {
  
  fit <- lm(
    logFC_5 ~ logFC_3,
    data = comparison[valid_logfc, ]
  )
  
  abline(fit, lwd = 2)
}

legend(
  "topleft",
  legend = c(
    "Zero logFC",
    "Identity line",
    "Linear fit"
  ),
  lty = c(2, 2, 1),
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 14. REPORT
# ------------------------------------------------------------

cat("\n13. GENERAZIONE REPORT\n")
cat("------------------------------------------------------------\n")

n_total <- nrow(comparison)
n_both <- nrow(significant_both)
n_concordant <- nrow(concordant)
n_discordant <- nrow(discordant)

pearson_all <- correlation_results$pearson_r[
  correlation_results$subset == "ALL_GENES"
]

spearman_all <- correlation_results$spearman_rho[
  correlation_results$subset == "ALL_GENES"
]

report <- c(
  "# TRANSCRIPTOME-X - Cross-Knockdown Concordance",
  "",
  "## Metodo",
  "",
  paste0("- Soglia FDR: ", FDR_THRESHOLD),
  "- Identificativo di integrazione: Ensembl stable ID",
  "- Confronti: siNRP1_3 vs siNT e siNRP1_5 vs siNT",
  "- Significatività valutata separatamente in ciascun confronto",
  "",
  "## Risultati",
  "",
  paste0("- Geni complessivi: ", n_total),
  paste0("- Significativi in entrambi: ", n_both),
  paste0("- Concordanti tra i due knockdown: ", n_concordant),
  paste0("- Discordanti tra i due knockdown: ", n_discordant),
  paste0("- Significativi solo in siNRP1_3: ", nrow(significant_3_only)),
  paste0("- Significativi solo in siNRP1_5: ", nrow(significant_5_only)),
  paste0("- Non significativi in entrambi: ", nrow(not_significant)),
  "",
  "## Correlazione dei logFC",
  "",
  paste0("- Pearson r: ", round(pearson_all, 4)),
  paste0("- Spearman rho: ", round(spearman_all, 4)),
  "",
  "## Interpretazione",
  "",
  "La concordanza indica che i due confronti presentano variazioni",
  "di espressione nella stessa direzione per i geni significativi",
  "in entrambi i knockdown.",
  "",
  "La discordanza indica direzioni opposte e richiede un approfondimento.",
  "La significatività in un solo confronto non dimostra necessariamente",
  "l'assenza di effetto nell'altro.",
  "",
  "La correlazione dei logFC misura l'associazione tra le stime",
  "di effetto e non costituisce da sola una prova di replicazione.",
  "",
  "I risultati sono osservazionali e non dimostrano causalità."
)

writeLines(
  report,
  file.path(
    output_dir,
    "TRANSCRIPTOME_X_CROSS_KNOCKDOWN_REPORT.md"
  )
)

# ------------------------------------------------------------
# 15. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n14. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "cross_knockdown_category_summary.csv",
  "significance_contingency_table.csv",
  "logFC_correlation_results.csv",
  "all_genes_cross_knockdown.csv",
  "genes_significant_both.csv",
  "genes_concordant_both_knockdowns.csv",
  "genes_discordant_both_knockdowns.csv",
  "genes_significant_siNRP1_3_only.csv",
  "genes_significant_siNRP1_5_only.csv",
  "genes_not_significant_both.csv",
  "concordant_gene_priority_ranking.csv",
  "R23_core_genes_cross_knockdown_validation.csv",
  "logFC_cross_knockdown_scatter.pdf",
  "TRANSCRIPTOME_X_CROSS_KNOCKDOWN_REPORT.md"
)

validation <- data.frame(
  file = output_files,
  exists = file.exists(
    file.path(output_dir, output_files)
  ),
  size_kb = round(
    file.info(
      file.path(output_dir, output_files)
    )$size / 1024,
    2
  )
)

print(validation)

if (!all(validation$exists)) {
  stop("Validazione output fallita.")
}

# ------------------------------------------------------------
# 16. RIEPILOGO FINALE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CROSS-KNOCKDOWN CONCORDANCE COMPLETED\n")
cat("============================================================\n")

cat("\nDISTRIBUZIONE RISULTATI\n")
cat("------------------------------------------------------------\n")

cat("Geni complessivi:", n_total, "\n")
cat("Significativi in entrambi:", n_both, "\n")
cat("Concordanti:", n_concordant, "\n")
cat("Discordanti:", n_discordant, "\n")
cat("Solo siNRP1_3:", nrow(significant_3_only), "\n")
cat("Solo siNRP1_5:", nrow(significant_5_only), "\n")
cat("Non significativi:", nrow(not_significant), "\n")

cat("\nCORRELAZIONE LOGFC\n")
cat("------------------------------------------------------------\n")

print(
  correlation_results,
  row.names = FALSE
)

cat("\nTOP 20 GENI CONCORDANTI\n")
cat("------------------------------------------------------------\n")

if (nrow(concordant) > 0) {
  
  print(
    head(
      concordant[
        ,
        c(
          "ensembl_id",
          "gene_name",
          "logFC_3",
          "FDR_3",
          "logFC_5",
          "FDR_5",
          "minimum_abs_logFC",
          "maximum_FDR"
        )
      ],
      20
    ),
    row.names = FALSE
  )
  
} else {
  
  cat("Nessun gene significativo e concordante in entrambi i confronti.\n")
}

cat(
  "\nCartella output:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nCROSS-KNOCKDOWN ANALYSIS VALIDATED.\n")
