
# ============================================================
# TRANSCRIPTOME-X
# STEP 6 - EXPRESSION HEATMAPS
# ============================================================
#
# Dataset: GSE266566
# Cell line: MDA-MB-231
# Comparisons:
#   1. siNRP1_3 vs siNT
#   2. siNRP1_5 vs siNT
#
# Method:
#   - TMM-normalized logCPM
#   - Top 50 DE genes per comparison (FDR < 0.05)
#   - Row-wise Z-score standardization
#   - Hierarchical clustering
#
# ============================================================


# 1. CONFIGURATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - EXPRESSION HEATMAPS\n")
cat("============================================================\n")

input_dge <- file.path(
  "data", "processed", "normalized",
  "GSE266566_edgeR_TMM.rds"
)

input_dir <- file.path(
  "data", "processed", "differential_expression",
  "annotated"
)

output_dir <- file.path(
  "data", "processed", "visualizations",
  "expression_heatmaps"
)

comparisons <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

n_genes <- 50

# Check working directory and input files

if (!file.exists(input_dge)) {
  stop("DGEList file not found: ", input_dge)
}

for (comparison in comparisons) {
  
  input_file <- file.path(
    input_dir,
    paste0(comparison, "_annotated.csv")
  )
  
  if (!file.exists(input_file)) {
    stop("Annotated results not found: ", input_file)
  }
}


# 2. PACKAGES
# ============================================================

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("Package edgeR is not installed.")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package pheatmap is not installed.")
}

library(edgeR)
library(pheatmap)


# 3. LOAD NORMALIZED DATA
# ============================================================

cat("\n1. IMPORTAZIONE DATI\n")
cat("------------------------------------------------------------\n")

dge <- readRDS(input_dge)

if (!inherits(dge, "DGEList")) {
  stop("The imported object is not an edgeR DGEList.")
}

if (ncol(dge$counts) != 9) {
  stop("Expected 9 samples. Found: ", ncol(dge$counts))
}

cat("Geni:", nrow(dge$counts), "\n")
cat("Campioni:", ncol(dge$counts), "\n")


# 4. SAMPLE ANNOTATION
# ============================================================

cat("\n2. ANNOTAZIONE CAMPIONI\n")
cat("------------------------------------------------------------\n")

sample_names <- colnames(dge$counts)

condition <- rep(NA_character_, length(sample_names))

condition[grepl("siNT", sample_names, fixed = TRUE)] <-
  "siNT"

condition[grepl("siNRP1.3", sample_names, fixed = TRUE)] <-
  "siNRP1_3"

condition[grepl("siNRP1.5", sample_names, fixed = TRUE)] <-
  "siNRP1_5"

if (anyNA(condition)) {
  stop(
    "Unable to identify all sample conditions. ",
    "Check sample names in the DGEList."
  )
}

replicate <- sub(
  ".*RNA\\.([0-9]+)$",
  "\\1",
  sample_names
)

if (any(replicate == sample_names)) {
  stop("Unable to extract replicate numbers from sample names.")
}

annotation_col <- data.frame(
  Condition = factor(
    condition,
    levels = c("siNT", "siNRP1_3", "siNRP1_5")
  ),
  Replicate = factor(replicate),
  row.names = sample_names
)

if (any(table(annotation_col$Condition) != 3)) {
  stop("Expected three biological replicates per condition.")
}

print(annotation_col)


# 5. CALCULATE NORMALIZED LOGCPM
# ============================================================

cat("\n3. CALCOLO LOGCPM NORMALIZZATI\n")
cat("------------------------------------------------------------\n")

logcpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

if (any(!is.finite(logcpm))) {
  stop("Non-finite values detected in normalized logCPM.")
}

cat("Dimensioni matrice logCPM:",
    nrow(logcpm), "x", ncol(logcpm), "\n")


# 6. OUTPUT DIRECTORY
# ============================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 7. FUNCTION TO GENERATE HEATMAP
# ============================================================

generate_heatmap <- function(comparison) {
  
  cat("\n")
  cat("============================================================\n")
  cat("CONFRONTO:", comparison, "\n")
  cat("============================================================\n")
  
  # Import annotated differential expression results
  
  input_file <- file.path(
    input_dir,
    paste0(comparison, "_annotated.csv")
  )
  
  results <- read.csv(
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
  
  if (!all(required_columns %in% names(results))) {
    stop("Missing required columns in: ", input_file)
  }
  
  # Select statistically significant genes
  
  significant <- results[
    !is.na(results$FDR) &
      results$FDR < 0.05,
    ,
    drop = FALSE
  ]
  
  significant <- significant[
    order(significant$FDR, significant$PValue),
    ,
    drop = FALSE
  ]
  
  if (nrow(significant) < n_genes) {
    stop(
      "Not enough significant genes for ",
      comparison
    )
  }
  
  selected <- head(significant, n_genes)
  
  cat("Geni significativi:", nrow(significant), "\n")
  cat("Geni selezionati:", nrow(selected), "\n")
  
  # Match selected genes with expression matrix
  
  selected_ids <- selected$gene_id
  
  if (anyDuplicated(selected_ids)) {
    stop("Duplicated gene IDs in differential expression results.")
  }
  
  if (!all(selected_ids %in% rownames(logcpm))) {
    stop("Some selected genes are missing from the logCPM matrix.")
  }
  
  expression_matrix <- logcpm[
    selected_ids,
    ,
    drop = FALSE
  ]
  
  # Remove genes with zero variance
  
  gene_sd <- apply(
    expression_matrix,
    1,
    sd
  )
  
  valid_genes <- is.finite(gene_sd) & gene_sd > 0
  
  expression_matrix <- expression_matrix[
    valid_genes,
    ,
    drop = FALSE
  ]
  
  selected <- selected[
    match(rownames(expression_matrix), selected$gene_id),
    ,
    drop = FALSE
  ]
  
  if (nrow(expression_matrix) == 0) {
    stop("No genes with non-zero expression variance.")
  }
  
  # Row-wise Z-score standardization
  
  zscore_matrix <- t(
    scale(
      t(expression_matrix)
    )
  )
  
  if (any(!is.finite(zscore_matrix))) {
    stop("Invalid values detected in Z-score matrix.")
  }
  
  # Gene labels
  
  gene_labels <- selected$gene_name
  
  missing_names <- is.na(gene_labels) |
    gene_labels == ""
  
  gene_labels[missing_names] <-
    selected$gene_id[missing_names]
  
  # Export heatmap data
  
  heatmap_data <- data.frame(
    gene_id = selected$gene_id,
    gene_name = gene_labels,
    zscore_matrix,
    check.names = FALSE
  )
  
  data_file <- file.path(
    output_dir,
    paste0(comparison, "_heatmap_data.csv")
  )
  
  write.csv(
    heatmap_data,
    data_file,
    row.names = FALSE
  )
  
  # Output paths
  
  png_file <- file.path(
    output_dir,
    paste0(comparison, "_heatmap.png")
  )
  
  pdf_file <- file.path(
    output_dir,
    paste0(comparison, "_heatmap.pdf")
  )
  
  # Heatmap title
  
  plot_title <- paste(
    "Expression Heatmap:",
    comparison,
    "\nTop",
    nrow(zscore_matrix),
    "DE genes"
  )
  
  # Generate PNG
  
  pheatmap(
    zscore_matrix,
    annotation_col = annotation_col,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_row = gene_labels,
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 9,
    border_color = NA,
    main = plot_title,
    filename = png_file,
    width = 11,
    height = 10
  )
  
  # Generate PDF
  
  pheatmap(
    zscore_matrix,
    annotation_col = annotation_col,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_row = gene_labels,
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 9,
    border_color = NA,
    main = plot_title,
    filename = pdf_file,
    width = 11,
    height = 10
  )
  
  cat("\nFile generati:\n")
  cat("PNG:", png_file, "\n")
  cat("PDF:", pdf_file, "\n")
  cat("CSV:", data_file, "\n")
  
  cat("\nVerifica NRP1:\n")
  
  nrp1 <- selected[
    selected$gene_name == "NRP1",
    ,
    drop = FALSE
  ]
  
  if (nrow(nrp1) > 0) {
    print(
      nrp1[, c("gene_id", "gene_name", "logFC", "FDR")]
    )
  } else {
    cat("NRP1 non rientra nei 50 geni selezionati.\n")
  }
  
  return(
    data.frame(
      comparison = comparison,
      significant_genes = nrow(significant),
      selected_genes = nrow(selected),
      genes_in_heatmap = nrow(zscore_matrix),
      nrp1_selected = nrow(nrp1) > 0,
      png_exists = file.exists(png_file),
      pdf_exists = file.exists(pdf_file),
      csv_exists = file.exists(data_file)
    )
  )
}


# 8. GENERATE BOTH HEATMAPS
# ============================================================

cat("\n4. GENERAZIONE HEATMAP\n")
cat("------------------------------------------------------------\n")

summary_list <- lapply(
  comparisons,
  generate_heatmap
)

heatmap_summary <- do.call(
  rbind,
  summary_list
)


# 9. EXPORT SUMMARY
# ============================================================

summary_file <- file.path(
  output_dir,
  "heatmap_summary.csv"
)

write.csv(
  heatmap_summary,
  summary_file,
  row.names = FALSE
)


# 10. FINAL VALIDATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RIEPILOGO EXPRESSION HEATMAPS\n")
cat("============================================================\n")

print(heatmap_summary)

if (!all(
  heatmap_summary$png_exists &
  heatmap_summary$pdf_exists &
  heatmap_summary$csv_exists
)) {
  stop("One or more output files were not generated.")
}

if (!file.exists(summary_file)) {
  stop("Heatmap summary file was not generated.")
}

cat("\nCONTROLLO FINALE SUPERATO\n")
cat("Tutti i file sono stati generati correttamente.\n")

cat("\nDirectory output:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
cat("STEP 6 COMPLETATO\n")
cat("============================================================\n")
