
# ============================================================
# TRANSCRIPTOME-X
# DIFFERENTIAL EXPRESSION - QUASI-LIKELIHOOD
# ============================================================

cat("\nTRANSCRIPTOME-X - DIFFERENTIAL EXPRESSION\n")
cat("==========================================\n")

# 1. PACCHETTI
# ------------------------------------------------------------

library(edgeR)
library(limma)

# 2. IMPORTAZIONE DGEList CON DISPERSIONI
# ------------------------------------------------------------

input_file <- file.path(
  "data",
  "processed",
  "differential_expression",
  "edgeR_DGEList_dispersion.rds"
)

if (!file.exists(input_file)) {
  stop("File DGEList non trovato: ", input_file)
}

dge <- readRDS(input_file)

cat("\n1. IMPORTAZIONE DATI\n")
cat("----------------------------\n")
cat("Geni:", nrow(dge), "\n")
cat("Campioni:", ncol(dge), "\n")

# 3. RICOSTRUZIONE DESIGN MATRIX
# ------------------------------------------------------------

group <- factor(
  dge$samples$group,
  levels = c("siNT", "siNRP1_3", "siNRP1_5")
)

design <- model.matrix(~ 0 + group)

colnames(design) <- c(
  "siNT",
  "siNRP1_3",
  "siNRP1_5"
)

rownames(design) <- colnames(dge)

cat("\n2. DESIGN MATRIX\n")
cat("----------------------------\n")
print(design)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix non valida.")
}

# 4. CONTRASTI
# ------------------------------------------------------------

contrasts_matrix <- makeContrasts(
  siNRP1_3_vs_siNT = siNRP1_3 - siNT,
  siNRP1_5_vs_siNT = siNRP1_5 - siNT,
  levels = design
)

cat("\n3. CONTRASTI\n")
cat("----------------------------\n")
print(contrasts_matrix)

# 5. STIMA DEL MODELLO QUASI-LIKELIHOOD
# ------------------------------------------------------------

cat("\n4. FIT DEL MODELLO QUASI-LIKELIHOOD\n")
cat("-----------------------------------\n")

fit <- glmQLFit(
  dge,
  design,
  robust = TRUE
)

cat("Modello stimato correttamente.\n")

# 6. SALVATAGGIO MODELLO
# ------------------------------------------------------------

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

saveRDS(
  fit,
  file.path(output_dir, "edgeR_QL_fit.rds")
)

cat("\nModello salvato.\n")

# 7. FINE FASE
# ------------------------------------------------------------

cat("\n==========================================\n")
cat("STEP 2 COMPLETATO\n")
cat("==========================================\n")

