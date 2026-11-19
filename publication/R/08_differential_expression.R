

# ============================================================
# TRANSCRIPTOME-X
# DIFFERENTIAL EXPRESSION ANALYSIS
# STEP 1 - DESIGN MATRIX AND DISPERSION ESTIMATION
# ============================================================

# 1. LIBRARIES
library(edgeR)
library(limma)

# 2. LOAD DATA
cat("\n1. IMPORTAZIONE DATI\n")
cat("------------------------------------------------------------\n")

dge <- readRDS(
  "data/processed/normalized/GSE266566_edgeR_TMM.rds"
)

cat("Geni:", nrow(dge), "\n")
cat("Campioni:", ncol(dge), "\n")

# 3. DEFINE EXPERIMENTAL GROUPS
cat("\n2. CONDIZIONI SPERIMENTALI\n")
cat("------------------------------------------------------------\n")

dge$samples$group <- factor(
  dge$samples$group,
  levels = c("siNT", "siNRP1_3", "siNRP1_5")
)

print(table(dge$samples$group))

# 4. DESIGN MATRIX
cat("\n3. DESIGN MATRIX\n")
cat("------------------------------------------------------------\n")

design <- model.matrix(
  ~ 0 + group,
  data = dge$samples
)

colnames(design) <- c(
  "siNT",
  "siNRP1_3",
  "siNRP1_5"
)

rownames(design) <- colnames(dge)

print(design)

# Validate design matrix
if (qr(design)$rank != ncol(design)) {
  stop("Design matrix non valida.")
}

cat("\nDesign matrix: VALIDATA\n")

# 5. DISPERSION ESTIMATION
cat("\n4. STIMA DELLE DISPERSIONI\n")
cat("------------------------------------------------------------\n")

dge <- estimateDisp(
  dge,
  design,
  robust = TRUE
)

cat(
  "Common dispersion:",
  dge$common.dispersion,
  "\n"
)

cat("\nRiepilogo tagwise dispersion:\n")

print(
  summary(dge$tagwise.dispersion)
)

# 6. BCV PLOT
cat("\n5. GENERAZIONE GRAFICO BCV\n")
cat("------------------------------------------------------------\n")

output_dir <- "data/processed/differential_expression"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(
  file.path(output_dir, "BCV_plot.pdf"),
  width = 8,
  height = 7
)

plotBCV(
  dge,
  main = "TRANSCRIPTOME-X - Biological Coefficient of Variation"
)

dev.off()

# 7. SAVE OBJECT
saveRDS(
  dge,
  file.path(
    output_dir,
    "edgeR_DGEList_dispersion.rds"
  )
)

cat("\nOggetto DGEList con dispersioni salvato.\n")

cat("\n============================================================\n")
cat("STEP 1 COMPLETATO\n")
cat("============================================================\n")

