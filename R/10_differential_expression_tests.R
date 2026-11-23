
# ============================================================
# TRANSCRIPTOME-X
# DIFFERENTIAL EXPRESSION TESTS
# edgeR - Quasi-Likelihood F-test
# ============================================================

cat("\nTRANSCRIPTOME-X - DIFFERENTIAL EXPRESSION TESTS\n")
cat("================================================\n")

# 1. PACCHETTI
# ------------------------------------------------------------

library(edgeR)

# 2. IMPORTAZIONE MODELLO
# ------------------------------------------------------------

input_file <- file.path(
  "data",
  "processed",
  "differential_expression",
  "edgeR_QL_fit.rds"
)

if (!file.exists(input_file)) {
  stop("Modello QL non trovato.")
}

fit <- readRDS(input_file)

output_dir <- file.path(
  "data",
  "processed",
  "differential_expression"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n1. MODELLO IMPORTATO\n")
cat("----------------------------\n")
cat("Geni:", nrow(fit), "\n")
cat("Campioni:", ncol(fit), "\n")

# 3. DEFINIZIONE CONTRASTI
# ------------------------------------------------------------

contrasts_matrix <- makeContrasts(
  siNRP1_3_vs_siNT = siNRP1_3 - siNT,
  siNRP1_5_vs_siNT = siNRP1_5 - siNT,
  levels = colnames(fit$design)
)

# 4. TEST 1: siNRP1_3 vs siNT
# ------------------------------------------------------------

cat("\n2. TEST siNRP1_3 vs siNT\n")
cat("----------------------------\n")

qlf_3 <- glmQLFTest(
  fit,
  contrast = contrasts_matrix[, "siNRP1_3_vs_siNT"]
)

results_3 <- topTags(
  qlf_3,
  n = Inf,
  sort.by = "none"
)$table

results_3$gene_id <- rownames(results_3)

results_3 <- results_3[
  ,
  c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
]

results_3 <- results_3[
  order(results_3$FDR),
]

write.csv(
  results_3,
  file.path(output_dir, "siNRP1_3_vs_siNT_all_genes.csv"),
  row.names = FALSE
)

cat("Geni testati:", nrow(results_3), "\n")
cat("FDR < 0.05:", sum(results_3$FDR < 0.05), "\n")
cat(
  "FDR < 0.05 e |logFC| >= 1:",
  sum(results_3$FDR < 0.05 & abs(results_3$logFC) >= 1),
  "\n"
)

# 5. TEST 2: siNRP1_5 vs siNT
# ------------------------------------------------------------

cat("\n3. TEST siNRP1_5 vs siNT\n")
cat("----------------------------\n")

qlf_5 <- glmQLFTest(
  fit,
  contrast = contrasts_matrix[, "siNRP1_5_vs_siNT"]
)

results_5 <- topTags(
  qlf_5,
  n = Inf,
  sort.by = "none"
)$table

results_5$gene_id <- rownames(results_5)

results_5 <- results_5[
  ,
  c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
]

results_5 <- results_5[
  order(results_5$FDR),
]

write.csv(
  results_5,
  file.path(output_dir, "siNRP1_5_vs_siNT_all_genes.csv"),
  row.names = FALSE
)

cat("Geni testati:", nrow(results_5), "\n")
cat("FDR < 0.05:", sum(results_5$FDR < 0.05), "\n")
cat(
  "FDR < 0.05 e |logFC| >= 1:",
  sum(results_5$FDR < 0.05 & abs(results_5$logFC) >= 1),
  "\n"
)

# 6. SALVATAGGIO RIEPILOGO
# ------------------------------------------------------------

summary_results <- data.frame(
  comparison = c(
    "siNRP1_3_vs_siNT",
    "siNRP1_5_vs_siNT"
  ),
  genes_tested = c(
    nrow(results_3),
    nrow(results_5)
  ),
  significant_FDR_005 = c(
    sum(results_3$FDR < 0.05),
    sum(results_5$FDR < 0.05)
  ),
  significant_FDR_005_logFC_1 = c(
    sum(results_3$FDR < 0.05 & abs(results_3$logFC) >= 1),
    sum(results_5$FDR < 0.05 & abs(results_5$logFC) >= 1)
  )
)

write.csv(
  summary_results,
  file.path(output_dir, "differential_expression_summary.csv"),
  row.names = FALSE
)

cat("\n4. RIEPILOGO\n")
cat("----------------------------\n")
print(summary_results)

# 7. SALVATAGGIO OGGETTI STATISTICI
# ------------------------------------------------------------

saveRDS(
  list(
    qlf_siNRP1_3 = qlf_3,
    qlf_siNRP1_5 = qlf_5,
    results_siNRP1_3 = results_3,
    results_siNRP1_5 = results_5
  ),
  file.path(output_dir, "edgeR_QLF_results.rds")
)

cat("\n================================================\n")
cat("STEP 3 COMPLETATO\n")
cat("Risultati salvati in:", output_dir, "\n")
cat("================================================\n")
