
# ============================================================
# TRANSCRIPTOME-X
# Script: 28_normalization_sensitivity.R
#
# NORMALIZATION SENSITIVITY ANALYSIS
#
# Metodi:
#   1. TMM
#   2. TMMwsp
#   3. RLE
#   4. upperquartile
#
# Modello:
#   edgeR quasi-likelihood
#
# Contrasti:
#   siNRP1_3 vs siNT
#   siNRP1_5 vs siNT
#
# Obiettivo:
# Valutare la robustezza dei risultati DE rispetto alla
# scelta del metodo di normalizzazione.
#
# Nota metodologica:
# - I conteggi RSEM originali non vengono arrotondati.
# - Il filtraggio viene eseguito una sola volta.
# - Il disegno sperimentale e' identico per tutti i metodi.
# - Non vengono ricostruiti effective lengths.
# - TMM e' il metodo di riferimento.
# ============================================================


# ============================================================
# 0. CONFIGURAZIONE
# ============================================================

project_dir <- "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"

counts_file <- file.path(
  project_dir,
  "data",
  "processed",
  "counts",
  "GSE266566_MDA_MB_231_selected_counts.tsv"
)

metadata_file <- file.path(
  project_dir,
  "data",
  "processed",
  "metadata",
  "sample_metadata.csv"
)

output_dir <- file.path(
  project_dir,
  "results",
  "normalization_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("NORMALIZATION SENSITIVITY ANALYSIS\n")
cat("============================================================\n")


# ============================================================
# 1. DIPENDENZE
# ============================================================

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR non e' installato.")
}

cat("\nedgeR version: ",
    as.character(packageVersion("edgeR")),
    "\n",
    sep = "")


# ============================================================
# 2. VERIFICA FILE INPUT
# ============================================================

if (!file.exists(counts_file)) {
  stop(
    "File dei conteggi non trovato:\n",
    counts_file
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "File metadata non trovato:\n",
    metadata_file
  )
}

cat("\nFile input verificati.\n")


# ============================================================
# 3. IMPORTAZIONE METADATA
# ============================================================

metadata <- read.csv(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample_id",
  "condition"
)

missing_columns <- setdiff(
  required_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti nel metadata: ",
    paste(missing_columns, collapse = ", ")
  )
}

metadata$sample_id <- as.character(metadata$sample_id)
metadata$condition <- as.character(metadata$condition)

if (anyNA(metadata$sample_id) ||
    any(metadata$sample_id == "")) {
  stop("Il metadata contiene sample_id mancanti o vuoti.")
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Il metadata contiene sample_id duplicati.")
}

if (anyNA(metadata$condition) ||
    any(metadata$condition == "")) {
  stop("Il metadata contiene condizioni mancanti o vuote.")
}

expected_groups <- c(
  "siNT",
  "siNRP1_3",
  "siNRP1_5"
)

if (!setequal(unique(metadata$condition), expected_groups)) {
  stop(
    "Le condizioni nel metadata non corrispondono a quelle attese.\n",
    "Trovate: ",
    paste(unique(metadata$condition), collapse = ", ")
  )
}

if (nrow(metadata) != 9) {
  stop(
    "Attesi 9 campioni, trovati ",
    nrow(metadata),
    "."
  )
}

if (!all(table(factor(
  metadata$condition,
  levels = expected_groups
)) == 3)) {
  stop("Attesi 3 replicati biologici per ciascuna condizione.")
}

metadata$condition <- factor(
  metadata$condition,
  levels = expected_groups
)

cat("\nDistribuzione dei campioni:\n")
print(table(metadata$condition))


# ============================================================
# 4. IMPORTAZIONE MATRICE DEI CONTEGGI
# ============================================================

counts_df <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(counts_df) < 4) {
  stop(
    "La matrice deve contenere almeno 3 colonne di annotazione ",
    "e una colonna di conteggi."
  )
}

# Identificazione Gene ID.
gene_id_column <- "Ensembl.114.Gene.ID"

if (!gene_id_column %in% colnames(counts_df)) {
  stop(
    "Colonna Gene ID non trovata: ",
    gene_id_column
  )
}

gene_id <- as.character(
  counts_df[[gene_id_column]]
)

if (anyNA(gene_id) || any(gene_id == "")) {
  stop("Sono presenti Gene ID mancanti o vuoti.")
}

if (anyDuplicated(gene_id)) {
  stop("Sono presenti Gene ID duplicati.")
}

# Le prime tre colonne contengono le annotazioni.
gene_info <- counts_df[, 1:3, drop = FALSE]

# Le colonne dei campioni devono corrispondere al metadata.
sample_ids <- metadata$sample_id

missing_samples <- setdiff(
  sample_ids,
  colnames(counts_df)
)

if (length(missing_samples) > 0) {
  stop(
    "Campioni mancanti nella matrice:\n",
    paste(missing_samples, collapse = "\n")
  )
}

# Estrazione esplicita dei nove campioni.
count_df <- counts_df[
  ,
  sample_ids,
  drop = FALSE
]

# Conversione numerica robusta.
count_list <- lapply(
  count_df,
  function(x) {
    suppressWarnings(
      as.numeric(as.character(x))
    )
  }
)

conversion_errors <- vapply(
  count_list,
  function(x) anyNA(x),
  logical(1)
)

if (any(conversion_errors)) {
  stop(
    "Conversione numerica non riuscita per le colonne: ",
    paste(names(conversion_errors)[conversion_errors],
          collapse = ", ")
  )
}

counts <- do.call(
  cbind,
  count_list
)

storage.mode(counts) <- "double"

rownames(counts) <- gene_id
colnames(counts) <- sample_ids

if (!identical(colnames(counts), metadata$sample_id)) {
  stop("L'ordine dei campioni non corrisponde al metadata.")
}

if (anyNA(counts)) {
  stop("La matrice contiene valori NA.")
}

if (any(!is.finite(counts))) {
  stop("La matrice contiene valori non finiti.")
}

if (any(counts < 0)) {
  stop("La matrice contiene conteggi negativi.")
}

if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("La matrice dei conteggi e' vuota.")
}

cat("\nMatrice dei conteggi:\n")
cat("Geni:", nrow(counts), "\n")
cat("Campioni:", ncol(counts), "\n")
cat(
  "Conteggi frazionari:",
  sum(counts != round(counts)),
  "\n"
)


# ============================================================
# 5. CONTROLLO IDENTITA' DEI CAMPIONI
# ============================================================

sample_audit <- data.frame(
  sample_id = metadata$sample_id,
  condition = as.character(metadata$condition),
  count_column = colnames(counts),
  total_counts = colSums(counts),
  stringsAsFactors = FALSE
)

if (!identical(
  sample_audit$sample_id,
  sample_audit$count_column
)) {
  stop("Fallimento del controllo di allineamento campioni.")
}

write.csv(
  sample_audit,
  file.path(output_dir, "sample_alignment_audit.csv"),
  row.names = FALSE
)


# ============================================================
# 6. CREAZIONE DGEList
# ============================================================

dge <- edgeR::DGEList(
  counts = counts,
  samples = metadata,
  genes = gene_info
)

dge$samples$group <- metadata$condition


# ============================================================
# 7. FILTRAGGIO
# ============================================================

cat("\nEsecuzione filterByExpr...\n")

keep <- edgeR::filterByExpr(
  dge,
  group = metadata$condition
)

if (length(keep) != nrow(counts)) {
  stop("Dimensione del vettore di filtraggio non valida.")
}

if (!any(keep)) {
  stop("Nessun gene mantenuto dal filtraggio.")
}

filtered_counts <- counts[
  keep,
  ,
  drop = FALSE
]

filtered_gene_info <- gene_info[
  keep,
  ,
  drop = FALSE
]

if (nrow(filtered_counts) < 2) {
  stop("Numero insufficiente di geni dopo il filtraggio.")
}

filter_summary <- data.frame(
  genes_initial = nrow(counts),
  genes_retained = nrow(filtered_counts),
  genes_removed = sum(!keep),
  samples = ncol(counts),
  stringsAsFactors = FALSE
)

write.csv(
  filter_summary,
  file.path(output_dir, "filtering_summary.csv"),
  row.names = FALSE
)

cat("\nFiltraggio completato.\n")
cat("Geni iniziali:", nrow(counts), "\n")
cat("Geni mantenuti:", nrow(filtered_counts), "\n")
cat("Geni esclusi:", sum(!keep), "\n")


# ============================================================
# 8. DISEGNO SPERIMENTALE
# ============================================================

group <- metadata$condition

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- expected_groups

rownames(design) <- metadata$sample_id

if (qr(design)$rank != ncol(design)) {
  stop("La matrice del disegno non e' a rango pieno.")
}

if (nrow(design) != ncol(filtered_counts)) {
  stop("Dimensioni incompatibili tra design e matrice dei conteggi.")
}

cat("\nDesign matrix:\n")
print(design)


# ============================================================
# 9. DEFINIZIONE CONTRASTI
# ============================================================

contrast_definitions <- list(
  
  siNRP1_3_vs_siNT = c(
    siNT = -1,
    siNRP1_3 = 1,
    siNRP1_5 = 0
  ),
  
  siNRP1_5_vs_siNT = c(
    siNT = -1,
    siNRP1_3 = 0,
    siNRP1_5 = 1
  )
)

contrasts_list <- lapply(
  contrast_definitions,
  function(x) {
    
    if (!setequal(names(x), colnames(design))) {
      stop("Nomi dei coefficienti del contrasto non validi.")
    }
    
    matrix(
      as.numeric(x[colnames(design)]),
      ncol = 1,
      dimnames = list(
        colnames(design),
        "contrast"
      )
    )
  }
)


# ============================================================
# 10. METODI DI NORMALIZZAZIONE
# ============================================================

normalization_methods <- c(
  "TMM",
  "TMMwsp",
  "RLE",
  "upperquartile"
)

reference_method <- "TMM"

results_by_method <- list()
normalization_factors <- list()
effective_library_sizes <- list()


# ============================================================
# 11. FUNZIONE PER CORRELAZIONE SPEARMAN
# ============================================================

safe_spearman <- function(x, y) {
  
  valid <- is.finite(x) & is.finite(y)
  
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  if (length(unique(x)) < 2 ||
      length(unique(y)) < 2) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )
}


# ============================================================
# 12. ANALISI DIFFERENZIALE
# ============================================================

for (method in normalization_methods) {
  
  cat("\n")
  cat("============================================================\n")
  cat("NORMALIZZAZIONE:", method, "\n")
  cat("============================================================\n")
  
  # DGEList indipendente per ciascun metodo.
  y <- edgeR::DGEList(
    counts = filtered_counts,
    samples = metadata,
    genes = filtered_gene_info
  )
  
  y$samples$group <- group
  
  # Normalizzazione.
  y <- edgeR::calcNormFactors(
    y,
    method = method
  )
  
  if (any(!is.finite(y$samples$norm.factors)) ||
      any(y$samples$norm.factors <= 0)) {
    stop(
      "Fattori di normalizzazione non validi per ",
      method
    )
  }
  
  current_effective_lib_size <-
    y$samples$lib.size * y$samples$norm.factors
  
  if (any(!is.finite(current_effective_lib_size)) ||
      any(current_effective_lib_size <= 0)) {
    stop(
      "Effective library size non valido per ",
      method
    )
  }
  
  normalization_factors[[method]] <- data.frame(
    sample_id = metadata$sample_id,
    condition = as.character(metadata$condition),
    method = method,
    library_size = y$samples$lib.size,
    norm_factor = y$samples$norm.factors,
    effective_library_size = current_effective_lib_size,
    stringsAsFactors = FALSE
  )
  
  effective_library_sizes[[method]] <-
    current_effective_lib_size
  
  cat("\nFattori di normalizzazione:\n")
  
  print(
    normalization_factors[[method]]
  )
  
  # Stima dispersione.
  y <- edgeR::estimateDisp(
    y,
    design,
    robust = TRUE
  )
  
  # Modello quasi-likelihood.
  fit <- edgeR::glmQLFit(
    y,
    design,
    robust = TRUE
  )
  
  method_results <- list()
  
  for (contrast_name in names(contrasts_list)) {
    
    cat("\nContrasto:", contrast_name, "\n")
    
    test <- edgeR::glmQLFTest(
      fit,
      contrast = contrasts_list[[contrast_name]]
    )
    
    result_table <- edgeR::topTags(
      test,
      n = Inf,
      sort.by = "none"
    )$table
    
    result_table$gene_id <- rownames(result_table)
    
    result_table <- result_table[
      ,
      c(
        "gene_id",
        setdiff(
          colnames(result_table),
          "gene_id"
        )
      ),
      drop = FALSE
    ]
    
    # Controllo risultati.
    if (nrow(result_table) != nrow(filtered_counts)) {
      stop(
        "Numero di risultati inatteso per ",
        method,
        " - ",
        contrast_name
      )
    }
    
    method_results[[contrast_name]] <- result_table
    
    output_file <- file.path(
      output_dir,
      paste0(
        "DE_",
        method,
        "_",
        contrast_name,
        ".csv"
      )
    )
    
    write.csv(
      result_table,
      output_file,
      row.names = FALSE
    )
    
    significant_count <- sum(
      !is.na(result_table$FDR) &
        result_table$FDR < 0.05
    )
    
    cat("Geni testati:", nrow(result_table), "\n")
    cat(
      "Geni significativi FDR < 0.05:",
      significant_count,
      "\n"
    )
  }
  
  results_by_method[[method]] <- method_results
  
  # Salvataggio oggetti del modello.
  saveRDS(
    list(
      method = method,
      dge = y,
      fit = fit,
      results = method_results
    ),
    file.path(
      output_dir,
      paste0("model_", method, ".rds")
    )
  )
  
  cat("\nMetodo completato:", method, "\n")
}


# ============================================================
# 13. ESPORTAZIONE FATTORI DI NORMALIZZAZIONE
# ============================================================

normalization_factor_table <- do.call(
  rbind,
  normalization_factors
)

rownames(normalization_factor_table) <- NULL

write.csv(
  normalization_factor_table,
  file.path(
    output_dir,
    "normalization_factors_all_methods.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. CONFRONTO FATTORI TRA METODI
# ============================================================

factor_wide <- data.frame(
  sample_id = metadata$sample_id,
  stringsAsFactors = FALSE
)

for (method in normalization_methods) {
  
  current <- normalization_factors[[method]]
  
  factor_wide[[paste0("norm_factor_", method)]] <-
    current$norm_factor
  
  factor_wide[[paste0("effective_library_size_", method)]] <-
    current$effective_library_size
}

write.csv(
  factor_wide,
  file.path(
    output_dir,
    "normalization_factors_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. CONFRONTO DIFFERENTIAL EXPRESSION CON TMM
# ============================================================

comparison_results <- list()

for (contrast_name in names(contrasts_list)) {
  
  reference_results <- results_by_method[[
    reference_method
  ]][[contrast_name]]
  
  reference_sig <- reference_results$gene_id[
    !is.na(reference_results$FDR) &
      reference_results$FDR < 0.05
  ]
  
  for (method in setdiff(
    normalization_methods,
    reference_method
  )) {
    
    alternative_results <- results_by_method[[
      method
    ]][[contrast_name]]
    
    alternative_sig <- alternative_results$gene_id[
      !is.na(alternative_results$FDR) &
        alternative_results$FDR < 0.05
    ]
    
    common <- merge(
      reference_results[
        ,
        c("gene_id", "logFC", "FDR"),
        drop = FALSE
      ],
      alternative_results[
        ,
        c("gene_id", "logFC", "FDR"),
        drop = FALSE
      ],
      by = "gene_id",
      suffixes = c("_TMM", "_alternative")
    )
    
    overlap <- intersect(
      reference_sig,
      alternative_sig
    )
    
    union_genes <- union(
      reference_sig,
      alternative_sig
    )
    
    jaccard <- if (length(union_genes) == 0) {
      NA_real_
    } else {
      length(overlap) / length(union_genes)
    }
    
    spearman_correlation <- safe_spearman(
      common$logFC_TMM,
      common$logFC_alternative
    )
    
    if (length(overlap) == 0) {
      
      direction_concordance <- NA_real_
      
    } else {
      
      ref_index <- match(
        overlap,
        reference_results$gene_id
      )
      
      alt_index <- match(
        overlap,
        alternative_results$gene_id
      )
      
      ref_direction <- sign(
        reference_results$logFC[ref_index]
      )
      
      alt_direction <- sign(
        alternative_results$logFC[alt_index]
      )
      
      direction_concordance <- mean(
        ref_direction == alt_direction,
        na.rm = TRUE
      )
      
      if (!is.finite(direction_concordance)) {
        direction_concordance <- NA_real_
      }
    }
    
    comparison_results[[
      paste(contrast_name, method, sep = "__")
    ]] <- data.frame(
      contrast = contrast_name,
      reference = reference_method,
      alternative = method,
      genes_tested = nrow(common),
      significant_TMM = length(reference_sig),
      significant_alternative = length(alternative_sig),
      significant_overlap = length(overlap),
      Jaccard_index = jaccard,
      Spearman_logFC = spearman_correlation,
      direction_concordance_shared_significant =
        direction_concordance,
      stringsAsFactors = FALSE
    )
  }
}

comparison_table <- do.call(
  rbind,
  comparison_results
)

rownames(comparison_table) <- NULL

write.csv(
  comparison_table,
  file.path(
    output_dir,
    "normalization_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. CONFRONTO NUMERO GENI SIGNIFICATIVI
# ============================================================

significant_summary <- do.call(
  rbind,
  lapply(
    normalization_methods,
    function(method) {
      
      do.call(
        rbind,
        lapply(
          names(contrasts_list),
          function(contrast_name) {
            
            result <- results_by_method[[
              method
            ]][[contrast_name]]
            
            data.frame(
              method = method,
              contrast = contrast_name,
              genes_tested = nrow(result),
              significant_FDR_0_05 = sum(
                !is.na(result$FDR) &
                  result$FDR < 0.05
              ),
              significant_FDR_0_01 = sum(
                !is.na(result$FDR) &
                  result$FDR < 0.01
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

rownames(significant_summary) <- NULL

write.csv(
  significant_summary,
  file.path(
    output_dir,
    "significant_genes_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. SALVATAGGIO RISULTATI COMPLETI
# ============================================================

complete_results <- list(
  project = "TRANSCRIPTOME-X",
  analysis = "Normalization sensitivity",
  reference_method = reference_method,
  normalization_methods = normalization_methods,
  sample_metadata = metadata,
  sample_alignment = sample_audit,
  design = design,
  filtered_gene_ids = rownames(filtered_counts),
  filtering_summary = filter_summary,
  normalization_factors = normalization_factors,
  effective_library_sizes = effective_library_sizes,
  differential_expression = results_by_method,
  comparison = comparison_table,
  significant_summary = significant_summary
)

saveRDS(
  complete_results,
  file.path(
    output_dir,
    "normalization_sensitivity_complete.rds"
  )
)


# ============================================================
# 18. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 19. REPORT FINALE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("NORMALIZATION SENSITIVITY - FINAL REPORT\n")
cat("============================================================\n")

cat("\nCampioni analizzati:\n")
print(table(metadata$condition))

cat("\nGeni iniziali:", nrow(counts), "\n")
cat("Geni mantenuti:", nrow(filtered_counts), "\n")

cat("\nMetodi utilizzati:\n")
print(normalization_methods)

cat("\nNumero di geni significativi:\n")
print(significant_summary, row.names = FALSE)

cat("\nConfronto rispetto a TMM:\n")
print(comparison_table, row.names = FALSE)

cat("\nDirectory risultati:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat("\nFile generati:\n")

print(
  list.files(
    output_dir,
    full.names = FALSE
  )
)

cat("\n")
cat("ANALISI COMPLETATA.\n")
cat("============================================================\n")
