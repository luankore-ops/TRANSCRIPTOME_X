
# ============================================================
# TRANSCRIPTOME-X
# Script: 29_normalization_gene_stability.R
#
# GENE-LEVEL NORMALIZATION STABILITY AUDIT
#
# Input:
#   Risultati DE prodotti da 28_normalization_sensitivity.R
#
# Metodi:
#   TMM
#   TMMwsp
#   RLE
#   upperquartile
#
# Contrasti:
#   siNRP1_3_vs_siNT
#   siNRP1_5_vs_siNT
#
# Obiettivi:
#   1. Identificare i geni significativi con tutti i metodi.
#   2. Identificare i geni significativi solo con alcuni metodi.
#   3. Identificare i geni mai significativi.
#   4. Quantificare la variabilita' di logFC e FDR.
#   5. Verificare la concordanza delle direzioni.
#   6. Individuare i geni vicini alla soglia FDR.
#
# Soglia DE:
#   FDR < 0.05
#
# Soglie descrittive aggiuntive:
#   FDR tra 0.04 e 0.06: prossimita' alla soglia.
#   Range logFC >= 0.5: variabilita' esplorativa.
#
# ATTENZIONE:
# Le soglie aggiuntive sono descrittive, non costituiscono
# criteri biologici universali.
#
# Nessun modello viene ricalcolato.
# Nessun dato originale viene modificato.
# ============================================================


# ============================================================
# 0. CONFIGURAZIONE
# ============================================================

project_dir <- "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"

input_dir <- file.path(
  project_dir,
  "results",
  "normalization_sensitivity"
)

output_dir <- file.path(
  project_dir,
  "results",
  "normalization_gene_stability"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

methods <- c(
  "TMM",
  "TMMwsp",
  "RLE",
  "upperquartile"
)

contrasts_to_analyze <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

fdr_threshold <- 0.05

fdr_near_lower <- 0.04
fdr_near_upper <- 0.06

logfc_range_threshold <- 0.5


cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("NORMALIZATION GENE STABILITY AUDIT\n")
cat("============================================================\n")


# ============================================================
# 1. VERIFICA DIRECTORY
# ============================================================

if (!dir.exists(input_dir)) {
  stop(
    "Directory dei risultati non trovata:\n",
    input_dir,
    "\nEseguire prima 28_normalization_sensitivity.R."
  )
}

cat("\nDirectory input verificata.\n")


# ============================================================
# 2. VERIFICA FILE INPUT
# ============================================================

input_files <- expand.grid(
  method = methods,
  contrast = contrasts_to_analyze,
  stringsAsFactors = FALSE
)

input_files$filename <- paste0(
  "DE_",
  input_files$method,
  "_",
  input_files$contrast,
  ".csv"
)

input_files$path <- file.path(
  input_dir,
  input_files$filename
)

missing_files <- input_files$filename[
  !file.exists(input_files$path)
]

if (length(missing_files) > 0) {
  stop(
    "File DE mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat(
  "File DE verificati:",
  nrow(input_files),
  "\n"
)


# ============================================================
# 3. FUNZIONE IMPORTAZIONE E VALIDAZIONE
# ============================================================

read_de_result <- function(file_path, method_name) {
  
  x <- read.csv(
    file_path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "logFC",
    "FDR"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Colonne mancanti in ",
      basename(file_path),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  x$gene_id <- as.character(x$gene_id)
  
  x$logFC <- suppressWarnings(
    as.numeric(as.character(x$logFC))
  )
  
  x$FDR <- suppressWarnings(
    as.numeric(as.character(x$FDR))
  )
  
  if (anyNA(x$gene_id) || any(x$gene_id == "")) {
    stop(
      "Gene ID mancanti in ",
      basename(file_path)
    )
  }
  
  if (anyDuplicated(x$gene_id)) {
    stop(
      "Gene ID duplicati in ",
      basename(file_path)
    )
  }
  
  if (any(!is.finite(x$logFC))) {
    stop(
      "Valori logFC non validi in ",
      basename(file_path)
    )
  }
  
  invalid_fdr <- !is.na(x$FDR) & (
    !is.finite(x$FDR) |
      x$FDR < 0 |
      x$FDR > 1
  )
  
  if (any(invalid_fdr)) {
    stop(
      "Valori FDR non validi in ",
      basename(file_path)
    )
  }
  
  if (all(is.na(x$FDR))) {
    stop(
      "Tutti i valori FDR sono mancanti in ",
      basename(file_path)
    )
  }
  
  result <- x[
    ,
    c("gene_id", "logFC", "FDR"),
    drop = FALSE
  ]
  
  names(result)[2:3] <- c(
    paste0("logFC_", method_name),
    paste0("FDR_", method_name)
  )
  
  result
}


# ============================================================
# 4. FUNZIONE RANGE ROBUSTO
# ============================================================

safe_range <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  diff(range(x))
}


# ============================================================
# 5. FUNZIONE CORRELAZIONE SPEARMAN
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
# 6. ANALISI DI UN CONTRASTO
# ============================================================

analyze_contrast <- function(contrast_name) {
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("CONTRASTO:", contrast_name, "\n")
  cat("------------------------------------------------------------\n")
  
  # ----------------------------------------------------------
  # 6.1 Importazione risultati
  # ----------------------------------------------------------
  
  result_list <- list()
  
  for (method in methods) {
    
    filename <- paste0(
      "DE_",
      method,
      "_",
      contrast_name,
      ".csv"
    )
    
    filepath <- file.path(
      input_dir,
      filename
    )
    
    result_list[[method]] <- read_de_result(
      filepath,
      method
    )
    
    cat(
      method,
      ":",
      nrow(result_list[[method]]),
      "geni\n"
    )
  }
  
  
  # ----------------------------------------------------------
  # 6.2 Verifica identita' dei geni
  # ----------------------------------------------------------
  
  reference_gene_ids <- result_list[["TMM"]]$gene_id
  
  for (method in methods) {
    
    current_gene_ids <- result_list[[method]]$gene_id
    
    if (!setequal(
      reference_gene_ids,
      current_gene_ids
    )) {
      stop(
        "Gli insiemi di Gene ID non coincidono tra TMM e ",
        method,
        " nel contrasto ",
        contrast_name
      )
    }
  }
  
  if (length(reference_gene_ids) !=
      nrow(result_list[["TMM"]])) {
    stop("Insieme di Gene ID TMM non valido.")
  }
  
  cat(
    "\nGene ID comuni a tutti i metodi:",
    length(reference_gene_ids),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 6.3 Costruzione matrice comparativa
  # ----------------------------------------------------------
  
  stability <- data.frame(
    gene_id = reference_gene_ids,
    stringsAsFactors = FALSE
  )
  
  for (method in methods) {
    
    current <- result_list[[method]]
    
    index <- match(
      reference_gene_ids,
      current$gene_id
    )
    
    if (anyNA(index)) {
      stop(
        "Errore nell'allineamento dei geni per ",
        method
      )
    }
    
    stability[[paste0("logFC_", method)]] <-
      current$logFC[index]
    
    stability[[paste0("FDR_", method)]] <-
      current$FDR[index]
  }
  
  
  # ----------------------------------------------------------
  # 6.4 Matrici logFC e FDR
  # ----------------------------------------------------------
  
  logfc_matrix <- as.matrix(
    stability[
      ,
      paste0("logFC_", methods),
      drop = FALSE
    ]
  )
  
  fdr_matrix <- as.matrix(
    stability[
      ,
      paste0("FDR_", methods),
      drop = FALSE
    ]
  )
  
  storage.mode(logfc_matrix) <- "double"
  storage.mode(fdr_matrix) <- "double"
  
  
  # ----------------------------------------------------------
  # 6.5 Classificazione significativita'
  # ----------------------------------------------------------
  
  significance_matrix <- !is.na(fdr_matrix) &
    fdr_matrix < fdr_threshold
  
  colnames(significance_matrix) <- methods
  
  stability$significant_methods_count <-
    rowSums(significance_matrix)
  
  stability$significant_methods <- apply(
    significance_matrix,
    1,
    function(x) {
      
      selected <- methods[x]
      
      if (length(selected) == 0) {
        return("NONE")
      }
      
      paste(
        selected,
        collapse = ";"
      )
    }
  )
  
  stability$stability_class <- ifelse(
    stability$significant_methods_count == length(methods),
    "SIGNIFICANT_ALL_METHODS",
    ifelse(
      stability$significant_methods_count == 0,
      "NOT_SIGNIFICANT_ANY_METHOD",
      "SIGNIFICANT_SUBSET_METHODS"
    )
  )
  
  
  # ----------------------------------------------------------
  # 6.6 Metriche logFC
  # ----------------------------------------------------------
  
  stability$logFC_min <- apply(
    logfc_matrix,
    1,
    min
  )
  
  stability$logFC_max <- apply(
    logfc_matrix,
    1,
    max
  )
  
  stability$logFC_range <- apply(
    logfc_matrix,
    1,
    safe_range
  )
  
  stability$max_abs_logFC_difference_vs_TMM <- apply(
    abs(
      sweep(
        logfc_matrix[
          ,
          methods != "TMM",
          drop = FALSE
        ],
        1,
        stability$logFC_TMM,
        FUN = "-"
      )
    ),
    1,
    max
  )
  
  # Concordanza della direzione tra tutti i metodi.
  stability$logFC_direction_concordant_all <- apply(
    logfc_matrix,
    1,
    function(x) {
      length(unique(sign(x))) == 1
    }
  )
  
  stability$logFC_range_flag <- (
    stability$logFC_range >= logfc_range_threshold
  )
  
  
  # ----------------------------------------------------------
  # 6.7 Metriche FDR
  # ----------------------------------------------------------
  
  stability$FDR_min <- apply(
    fdr_matrix,
    1,
    function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      min(x, na.rm = TRUE)
    }
  )
  
  stability$FDR_max <- apply(
    fdr_matrix,
    1,
    function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      max(x, na.rm = TRUE)
    }
  )
  
  stability$FDR_range <- apply(
    fdr_matrix,
    1,
    safe_range
  )
  
  stability$FDR_near_threshold <- apply(
    fdr_matrix,
    1,
    function(x) {
      any(
        !is.na(x) &
          x >= fdr_near_lower &
          x <= fdr_near_upper
      )
    }
  )
  
  stability$FDR_threshold_sensitive <- (
    stability$significant_methods_count > 0 &
      stability$significant_methods_count < length(methods)
  )
  
  
  # ----------------------------------------------------------
  # 6.8 Correlazioni logFC rispetto a TMM
  # ----------------------------------------------------------
  
  stability$Spearman_TMM_vs_TMMwsp <- NA_real_
  stability$Spearman_TMM_vs_RLE <- NA_real_
  stability$Spearman_TMM_vs_upperquartile <- NA_real_
  
  for (method in setdiff(methods, "TMM")) {
    
    correlation <- safe_spearman(
      stability$logFC_TMM,
      stability[[paste0("logFC_", method)]]
    )
    
    stability[[
      paste0("Spearman_TMM_vs_", method)
    ]] <- correlation
  }
  
  
  # ----------------------------------------------------------
  # 6.9 Riordino colonne
  # ----------------------------------------------------------
  
  first_columns <- c(
    "gene_id",
    "stability_class",
    "significant_methods_count",
    "significant_methods"
  )
  
  metric_columns <- c(
    "logFC_min",
    "logFC_max",
    "logFC_range",
    "max_abs_logFC_difference_vs_TMM",
    "logFC_direction_concordant_all",
    "logFC_range_flag",
    "FDR_min",
    "FDR_max",
    "FDR_range",
    "FDR_near_threshold",
    "FDR_threshold_sensitive"
  )
  
  method_columns <- unlist(
    lapply(
      methods,
      function(m) {
        c(
          paste0("logFC_", m),
          paste0("FDR_", m)
        )
      }
    )
  )
  
  correlation_columns <- paste0(
    "Spearman_TMM_vs_",
    setdiff(methods, "TMM")
  )
  
  stability <- stability[
    ,
    c(
      first_columns,
      metric_columns,
      method_columns,
      correlation_columns
    ),
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # 6.10 Esportazione tabella completa
  # ----------------------------------------------------------
  
  write.csv(
    stability,
    file.path(
      output_dir,
      paste0(
        "gene_stability_",
        contrast_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.11 Esportazione categorie
  # ----------------------------------------------------------
  
  categories <- c(
    "SIGNIFICANT_ALL_METHODS",
    "SIGNIFICANT_SUBSET_METHODS",
    "NOT_SIGNIFICANT_ANY_METHOD"
  )
  
  for (category in categories) {
    
    subset_data <- stability[
      stability$stability_class == category,
      ,
      drop = FALSE
    ]
    
    filename <- paste0(
      "genes_",
      category,
      "_",
      contrast_name,
      ".csv"
    )
    
    write.csv(
      subset_data,
      file.path(output_dir, filename),
      row.names = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # 6.12 Esportazione geni con variabilita' logFC
  # ----------------------------------------------------------
  
  variable_logfc <- stability[
    stability$logFC_range_flag,
    ,
    drop = FALSE
  ]
  
  variable_logfc <- variable_logfc[
    order(
      variable_logfc$logFC_range,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  write.csv(
    variable_logfc,
    file.path(
      output_dir,
      paste0(
        "genes_variable_logFC_",
        contrast_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.13 Esportazione geni vicini alla soglia FDR
  # ----------------------------------------------------------
  
  near_threshold <- stability[
    stability$FDR_near_threshold,
    ,
    drop = FALSE
  ]
  
  near_threshold <- near_threshold[
    order(
      near_threshold$FDR_min,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  write.csv(
    near_threshold,
    file.path(
      output_dir,
      paste0(
        "genes_near_FDR_threshold_",
        contrast_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.14 Correlazioni globali
  # ----------------------------------------------------------
  
  correlations <- data.frame(
    contrast = contrast_name,
    reference = "TMM",
    alternative = setdiff(methods, "TMM"),
    Spearman_logFC = vapply(
      setdiff(methods, "TMM"),
      function(method) {
        safe_spearman(
          stability$logFC_TMM,
          stability[[paste0("logFC_", method)]]
        )
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.15 Riepilogo categorie
  # ----------------------------------------------------------
  
  category_summary <- as.data.frame(
    table(
      factor(
        stability$stability_class,
        levels = categories
      )
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(category_summary) <- c(
    "stability_class",
    "gene_count"
  )
  
  category_summary$contrast <- contrast_name
  
  category_summary$percentage <- (
    category_summary$gene_count /
      nrow(stability)
  ) * 100
  
  category_summary <- category_summary[
    ,
    c(
      "contrast",
      "stability_class",
      "gene_count",
      "percentage"
    ),
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # 6.16 Riepilogo metriche
  # ----------------------------------------------------------
  
  metric_summary <- data.frame(
    contrast = contrast_name,
    genes_analyzed = nrow(stability),
    
    significant_all_methods = sum(
      stability$stability_class ==
        "SIGNIFICANT_ALL_METHODS"
    ),
    
    significant_subset_methods = sum(
      stability$stability_class ==
        "SIGNIFICANT_SUBSET_METHODS"
    ),
    
    not_significant_any_method = sum(
      stability$stability_class ==
        "NOT_SIGNIFICANT_ANY_METHOD"
    ),
    
    genes_near_FDR_threshold = sum(
      stability$FDR_near_threshold
    ),
    
    genes_logFC_range_ge_0_5 = sum(
      stability$logFC_range_flag
    ),
    
    genes_direction_concordant_all = sum(
      stability$logFC_direction_concordant_all
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 6.17 Report contrasto
  # ----------------------------------------------------------
  
  cat("\nClassificazione geni:\n")
  print(category_summary, row.names = FALSE)
  
  cat("\nMetriche riassuntive:\n")
  print(metric_summary, row.names = FALSE)
  
  cat("\nCorrelazioni logFC:\n")
  print(correlations, row.names = FALSE)
  
  cat(
    "\nGeni significativi con tutti i metodi:",
    metric_summary$significant_all_methods,
    "\n"
  )
  
  cat(
    "Geni significativi solo con alcuni metodi:",
    metric_summary$significant_subset_methods,
    "\n"
  )
  
  cat(
    "Geni mai significativi:",
    metric_summary$not_significant_any_method,
    "\n"
  )
  
  list(
    stability = stability,
    category_summary = category_summary,
    metric_summary = metric_summary,
    correlations = correlations
  )
}


# ============================================================
# 7. ESECUZIONE DEI DUE CONTRASTI
# ============================================================

analysis_results <- list()

for (contrast_name in contrasts_to_analyze) {
  
  analysis_results[[contrast_name]] <- analyze_contrast(
    contrast_name
  )
}


# ============================================================
# 8. UNIONE RIEPILOGHI
# ============================================================

category_summary_all <- do.call(
  rbind,
  lapply(
    analysis_results,
    function(x) x$category_summary
  )
)

metric_summary_all <- do.call(
  rbind,
  lapply(
    analysis_results,
    function(x) x$metric_summary
  )
)

correlation_summary_all <- do.call(
  rbind,
  lapply(
    analysis_results,
    function(x) x$correlations
  )
)

rownames(category_summary_all) <- NULL
rownames(metric_summary_all) <- NULL
rownames(correlation_summary_all) <- NULL


# ============================================================
# 9. ESPORTAZIONE RIEPILOGHI GLOBALI
# ============================================================

write.csv(
  category_summary_all,
  file.path(
    output_dir,
    "gene_stability_category_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  metric_summary_all,
  file.path(
    output_dir,
    "gene_stability_metric_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  correlation_summary_all,
  file.path(
    output_dir,
    "gene_stability_correlation_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. SALVATAGGIO OGGETTO COMPLETO
# ============================================================

complete_stability_results <- list(
  project = "TRANSCRIPTOME-X",
  analysis = "Normalization gene stability",
  methods = methods,
  contrasts = contrasts_to_analyze,
  fdr_threshold = fdr_threshold,
  fdr_near_threshold = c(
    fdr_near_lower,
    fdr_near_upper
  ),
  logfc_range_threshold = logfc_range_threshold,
  results = analysis_results,
  category_summary = category_summary_all,
  metric_summary = metric_summary_all,
  correlation_summary = correlation_summary_all
)

saveRDS(
  complete_stability_results,
  file.path(
    output_dir,
    "normalization_gene_stability_complete.rds"
  )
)


# ============================================================
# 11. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


# ============================================================
# 12. REPORT FINALE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("NORMALIZATION GENE STABILITY - FINAL REPORT\n")
cat("============================================================\n")

cat("\nRiepilogo classificazione:\n")
print(
  category_summary_all,
  row.names = FALSE
)

cat("\nRiepilogo metriche:\n")
print(
  metric_summary_all,
  row.names = FALSE
)

cat("\nCorrelazioni logFC:\n")
print(
  correlation_summary_all,
  row.names = FALSE
)

cat("\nDirectory output:\n")
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
