suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
})

cat("============================================================\n")
cat("TRANSCRIPTOME-X - EXPLORATORY QUALITY CONTROL\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS AND INPUT
# ------------------------------------------------------------

input_file <- "data/processed/normalized/GSE266566_edgeR_TMM.rds"
output_dir <- "data/processed/qc"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
    stop("Oggetto edgeR non trovato: ", input_file)
}

dge <- readRDS(input_file)

if (!inherits(dge, "DGEList")) {
    stop("Il file RDS non contiene un oggetto DGEList.")
}

# ------------------------------------------------------------
# 2. LOG-CPM
# ------------------------------------------------------------

cat("1. PREPARAZIONE DATI QC\n")
cat("------------------------------------------------------------\n")

logCPM <- cpm(
    dge,
    log = TRUE,
    prior.count = 2,
    normalized.lib.sizes = TRUE
)

if (any(!is.finite(logCPM))) {
    stop("Valori non finiti nella matrice logCPM.")
}

cat("Geni:", nrow(logCPM), "\n")
cat("Campioni:", ncol(logCPM), "\n")

# ------------------------------------------------------------
# 3. SAMPLE METADATA
# ------------------------------------------------------------

sample_metadata <- dge$samples

sample_metadata$sample_id <- rownames(sample_metadata)
sample_metadata$condition <- factor(
    as.character(sample_metadata$group),
    levels = c("siNT", "siNRP1_3", "siNRP1_5")
)

sample_metadata <- sample_metadata[
    match(colnames(logCPM), sample_metadata$sample_id),
    ,
    drop = FALSE
]

if (!identical(sample_metadata$sample_id, colnames(logCPM))) {
    stop("Ordine dei campioni non coerente.")
}

# ------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------

cat("\n2. PRINCIPAL COMPONENT ANALYSIS\n")
cat("------------------------------------------------------------\n")

pca <- prcomp(
    t(logCPM),
    center = TRUE,
    scale. = FALSE
)

variance <- pca$sdev^2 / sum(pca$sdev^2) * 100

pca_scores <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = sample_metadata$condition
)

cat("Varianza PC1:", round(variance[1], 2), "%\n")
cat("Varianza PC2:", round(variance[2], 2), "%\n")
cat(
    "Varianza cumulativa PC1 + PC2:",
    round(sum(variance[1:2]), 2),
    "%\n"
)

write.csv(
    pca_scores,
    file.path(output_dir, "PCA_scores.csv"),
    row.names = FALSE
)

write.csv(
    data.frame(
        component = paste0("PC", seq_along(variance)),
        variance_percent = variance
    ),
    file.path(output_dir, "PCA_variance.csv"),
    row.names = FALSE
)

p_pca <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = condition)
) +
    geom_point(size = 4) +
    geom_text(
        aes(label = sample_id),
        size = 2.7,
        vjust = -1,
        show.legend = FALSE
    ) +
    labs(
        title = "TRANSCRIPTOME-X - PCA",
        subtitle = "MDA-MB-231 | RSEM expected counts | TMM log-CPM",
        x = paste0("PC1 (", round(variance[1], 2), "%)"),
        y = paste0("PC2 (", round(variance[2], 2), "%)"),
        color = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
    )

ggsave(
    file.path(output_dir, "PCA_samples.png"),
    p_pca,
    width = 11,
    height = 7,
    dpi = 300
)

ggsave(
    file.path(output_dir, "PCA_samples.pdf"),
    p_pca,
    width = 11,
    height = 7
)

print(pca_scores)

# ------------------------------------------------------------
# 5. SAMPLE CORRELATION
# ------------------------------------------------------------

cat("\n3. SAMPLE CORRELATION\n")
cat("------------------------------------------------------------\n")

cor_matrix <- cor(
    logCPM,
    method = "pearson",
    use = "complete.obs"
)

write.csv(
    cor_matrix,
    file.path(output_dir, "sample_correlations.csv")
)

cat("Correlazione minima:",
    round(min(cor_matrix[upper.tri(cor_matrix)]), 4), "\n")

cat("Correlazione massima:",
    round(max(cor_matrix[upper.tri(cor_matrix)]), 4), "\n")

pheatmap(
    cor_matrix,
    main = "Sample Pearson Correlation",
    display_numbers = TRUE,
    number_format = "%.2f",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    filename = file.path(
        output_dir,
        "sample_correlation_heatmap.png"
    ),
    width = 10,
    height = 9
)

# ------------------------------------------------------------
# 6. SAMPLE DISTANCE
# ------------------------------------------------------------

cat("\n4. SAMPLE DISTANCE\n")
cat("------------------------------------------------------------\n")

distance_matrix <- as.matrix(
    dist(t(logCPM), method = "euclidean")
)

write.csv(
    distance_matrix,
    file.path(output_dir, "sample_distances.csv")
)

cat("Distanza minima:",
    round(min(distance_matrix[upper.tri(distance_matrix)]), 4), "\n")

cat("Distanza massima:",
    round(max(distance_matrix[upper.tri(distance_matrix)]), 4), "\n")

pheatmap(
    distance_matrix,
    main = "Euclidean Sample Distance",
    display_numbers = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    filename = file.path(
        output_dir,
        "sample_distance_heatmap.png"
    ),
    width = 10,
    height = 9
)

# ------------------------------------------------------------
# 7. FINAL VALIDATION
# ------------------------------------------------------------

cat("\n5. FINAL VALIDATION\n")
cat("------------------------------------------------------------\n")

expected_files <- c(
    "PCA_scores.csv",
    "PCA_variance.csv",
    "PCA_samples.png",
    "PCA_samples.pdf",
    "sample_correlations.csv",
    "sample_correlation_heatmap.png",
    "sample_distances.csv",
    "sample_distance_heatmap.png"
)

file_status <- file.exists(
    file.path(output_dir, expected_files)
)

for (i in seq_along(expected_files)) {
    cat(
        if (file_status[i]) "[OK] " else "[ERROR] ",
        expected_files[i],
        "\n"
    )
}

checks <- c(
    nrow(pca_scores) == 9,
    all(is.finite(pca_scores$PC1)),
    all(is.finite(pca_scores$PC2)),
    all(is.finite(cor_matrix)),
    all(is.finite(distance_matrix)),
    all(file_status)
)

if (!all(checks)) {
    stop("Validazione QC non superata.")
}

cat("\n============================================================\n")
cat("QC COMPLETATO\n")
cat("Campioni analizzati: 9\n")
cat("Geni analizzati:", nrow(logCPM), "\n")
cat("Campioni esclusi: 0\n")
cat("VALIDAZIONE: SUPERATA\n")
cat("============================================================\n")
