suppressPackageStartupMessages({
    library(edgeR)
})

cat("============================================================\n")
cat("TRANSCRIPTOME-X - FILTERING AND TMM NORMALIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

counts_file <- "data/processed/counts/GSE266566_MDA_MB_231_selected_counts.tsv"
metadata_file <- "data/processed/metadata/sample_metadata.csv"

output_dir <- "data/processed/normalized"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. IMPORT DATA
# ------------------------------------------------------------

cat("1. IMPORTAZIONE DATI\n")
cat("------------------------------------------------------------\n")

counts_data <- read.delim(
    counts_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

metadata <- read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sample_ids <- metadata$sample_id

if (!all(sample_ids %in% colnames(counts_data))) {
    stop("Campioni del metadata non corrispondenti alla matrice.")
}

counts <- as.matrix(
    counts_data[, sample_ids, drop = FALSE]
)

storage.mode(counts) <- "numeric"

gene_info <- counts_data[, 1:3]
rownames(gene_info) <- gene_info$Ensembl.114.Gene.ID

rownames(counts) <- gene_info$Ensembl.114.Gene.ID

metadata <- metadata[
    match(colnames(counts), metadata$sample_id),
    ,
    drop = FALSE
]

stopifnot(identical(colnames(counts), metadata$sample_id))

metadata$condition <- factor(
    metadata$condition,
    levels = c("siNT", "siNRP1_3", "siNRP1_5")
)

if (anyNA(metadata$condition)) {
    stop("Condizioni sperimentali non riconosciute.")
}

cat("Geni iniziali:", nrow(counts), "\n")
cat("Campioni:", ncol(counts), "\n")
print(table(metadata$condition))

# ------------------------------------------------------------
# 3. CREATE DGEList
# ------------------------------------------------------------

cat("\n2. CREAZIONE DGEList\n")
cat("------------------------------------------------------------\n")

dge <- DGEList(
    counts = counts,
    samples = metadata,
    genes = gene_info
)

dge$samples$group <- metadata$condition

# ------------------------------------------------------------
# 4. LOW-EXPRESSION FILTERING
# ------------------------------------------------------------

cat("\n3. FILTRAGGIO GENI POCO ESPRESSI\n")
cat("------------------------------------------------------------\n")

keep <- filterByExpr(
    dge,
    group = dge$samples$group
)

genes_before <- nrow(dge)
genes_after <- sum(keep)
genes_removed <- genes_before - genes_after

if (genes_after == 0) {
    stop("Nessun gene rimasto dopo il filtraggio.")
}

dge <- dge[
    keep,
    ,
    keep.lib.sizes = FALSE
]

cat("Geni iniziali:", genes_before, "\n")
cat("Geni conservati:", genes_after, "\n")
cat("Geni rimossi:", genes_removed, "\n")
cat(
    "Percentuale conservata:",
    round(100 * genes_after / genes_before, 2),
    "%\n"
)

# ------------------------------------------------------------
# 5. TMM NORMALIZATION
# ------------------------------------------------------------

cat("\n4. NORMALIZZAZIONE TMM\n")
cat("------------------------------------------------------------\n")

dge <- calcNormFactors(
    dge,
    method = "TMM"
)

dge$samples$effective_lib_size <-
    dge$samples$lib.size * dge$samples$norm.factors

print(
    dge$samples[
        ,
        c("group", "lib.size", "norm.factors", "effective_lib_size")
    ]
)

# ------------------------------------------------------------
# 6. VALIDATION
# ------------------------------------------------------------

cat("\n5. VALIDAZIONE\n")
cat("------------------------------------------------------------\n")

checks <- c(
    nrow(dge) > 0,
    ncol(dge) == 9,
    all(is.finite(dge$counts)),
    all(dge$counts >= 0),
    all(is.finite(dge$samples$norm.factors)),
    all(dge$samples$norm.factors > 0),
    all(dge$samples$effective_lib_size > 0),
    identical(colnames(dge), metadata$sample_id)
)

if (!all(checks)) {
    print(checks)
    stop("Validazione non superata.")
}

# ------------------------------------------------------------
# 7. SAVE RESULTS
# ------------------------------------------------------------

cat("\n6. SALVATAGGIO RISULTATI\n")
cat("------------------------------------------------------------\n")

saveRDS(
    dge,
    file.path(output_dir, "GSE266566_edgeR_TMM.rds")
)

write.csv(
    data.frame(
        sample_id = rownames(dge$samples),
        condition = dge$samples$group,
        library_size = dge$samples$lib.size,
        normalization_factor = dge$samples$norm.factors,
        effective_library_size = dge$samples$effective_lib_size
    ),
    file.path(output_dir, "TMM_normalization_factors.csv"),
    row.names = FALSE
)

write.csv(
    data.frame(
        metric = c(
            "genes_before_filtering",
            "genes_after_filtering",
            "genes_removed",
            "samples"
        ),
        value = c(
            genes_before,
            genes_after,
            genes_removed,
            ncol(dge)
        )
    ),
    file.path(output_dir, "filtering_summary.csv"),
    row.names = FALSE
)

# ------------------------------------------------------------
# 8. FINAL REPORT
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL REPORT\n")
cat("============================================================\n")

cat("Geni iniziali:", genes_before, "\n")
cat("Geni finali:", genes_after, "\n")
cat("Campioni:", ncol(dge), "\n")
cat("Metodo normalizzazione: TMM\n")
cat("Conteggi originali arrotondati: NO\n")
cat("Conteggi originali trasformati: NO\n")
cat("Oggetto DGEList salvato: SI\n")

cat("\nFile generati:\n")
cat("- GSE266566_edgeR_TMM.rds\n")
cat("- TMM_normalization_factors.csv\n")
cat("- filtering_summary.csv\n")

cat("\nVALIDAZIONE: SUPERATA\n")
cat("============================================================\n")
