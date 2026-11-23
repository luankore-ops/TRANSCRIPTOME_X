
# ============================================================
# TRANSCRIPTOME-X
# Script: 30_normalization_discordant_genes.R
#
# NORMALIZATION DISCORDANT GENES - COMPREHENSIVE AUDIT
#
# Obiettivi:
# 1. Identificare geni significativi solo con TMM.
# 2. Identificare geni significativi solo con metodi alternativi.
# 3. Identificare geni specifici di ciascun metodo.
# 4. Analizzare i profili di significativita' tra metodi.
# 5. Quantificare la distanza dalla soglia FDR = 0.05.
# 6. Confrontare FDR e logFC tra metodi.
# 7. Verificare la concordanza delle direzioni.
# 8. Annotare gli Ensembl ID con gene symbol, se disponibile.
# 9. Esportare tabelle complete e riepiloghi.
#
# Input:
# results/normalization_sensitivity/DE_<METHOD>_<CONTRAST>.csv
#
# Output:
# results/normalization_discordant_genes/
#
# Non ricalcola normalizzazione o modelli statistici.
# ============================================================


# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

project_dir <- "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"

input_dir <- file.path(
  project_dir,
  "results",
  "normalization_sensitivity"
)

output_dir <- file.path(
  project_dir,
  "results",
  "normalization_discordant_genes"
)

methods <- c(
  "TMM",
  "TMMwsp",
  "RLE",
  "upperquartile"
)

contrasts <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

fdr_threshold <- 0.05
near_fdr_window <- 0.01

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. FUNZIONI
# ------------------------------------------------------------

read_de_result <- function(method, contrast) {
  
  filename <- paste0(
    "DE_",
    method,
    "_",
    contrast,
    ".csv"
  )
  
  filepath <- file.path(
    input_dir,
    filename
  )
  
  if (!file.exists(filepath)) {
    stop(
      "File mancante: ",
      filepath
    )
  }
  
  data <- read.csv(
    filepath,
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
    colnames(data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Colonne mancanti in ",
      filename,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data$gene_id <- as.character(data$gene_id)
  
  data$logFC <- suppressWarnings(
    as.numeric(as.character(data$logFC))
  )
  
  data$FDR <- suppressWarnings(
    as.numeric(as.character(data$FDR))
  )
  
  if (anyNA(data$gene_id) ||
      any(data$gene_id == "")) {
    stop(
      "Gene ID mancanti o vuoti in ",
      filename
    )
  }
  
  if (anyDuplicated(data$gene_id) > 0) {
    stop(
      "Gene ID duplicati in ",
      filename
    )
  }
  
  if (any(!is.finite(data$logFC))) {
    stop(
      "Valori logFC non finiti in ",
      filename
    )
  }
  
  invalid_fdr <- !is.na(data$FDR) &
    (
      !is.finite(data$FDR) |
        data$FDR < 0 |
        data$FDR > 1
    )
  
  if (any(invalid_fdr)) {
    stop(
      "Valori FDR non validi in ",
      filename
    )
  }
  
  data[, c("gene_id", "logFC", "FDR")]
}


safe_min <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(x)
}


safe_max <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(x)
}


safe_range <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(x) - min(x)
}


safe_cor <- function(x, y) {
  
  valid <- is.finite(x) & is.finite(y)
  
  if (sum(valid) < 3) {
    return(NA_real_)
  }
  
  if (length(unique(x[valid])) < 2 ||
      length(unique(y[valid])) < 2) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(
      x[valid],
      y[valid],
      method = "spearman"
    )
  )
}


# ------------------------------------------------------------
# 3. CONTROLLO INPUT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("NORMALIZATION DISCORDANT GENES AUDIT\n")
cat("============================================================\n\n")

if (!dir.exists(input_dir)) {
  stop(
    "Directory input non trovata: ",
    input_dir
  )
}

expected_files <- expand.grid(
  method = methods,
  contrast = contrasts,
  stringsAsFactors = FALSE
)

expected_files$filename <- paste0(
  "DE_",
  expected_files$method,
  "_",
  expected_files$contrast,
  ".csv"
)

missing_files <- expected_files$filename[
  !file.exists(
    file.path(
      input_dir,
      expected_files$filename
    )
  )
]

if (length(missing_files) > 0) {
  stop(
    "File DE mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Directory input verificata.\n")
cat("File DE disponibili:", nrow(expected_files), "\n")
cat("Metodi:", paste(methods, collapse = ", "), "\n")
cat("Contrasti:", paste(contrasts, collapse = ", "), "\n")


# ------------------------------------------------------------
# 4. CARICAMENTO ANNOTAZIONE
# ------------------------------------------------------------

annotation_available <- FALSE
annotation_table <- NULL

if (
  requireNamespace("AnnotationDbi", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
) {
  
  annotation_available <- TRUE
  
  cat("\nAnnotazione org.Hs.eg.db disponibile.\n")
  
} else {
  
  cat(
    "\nAnnotazione org.Hs.eg.db non disponibile:",
    "verra' mantenuto il gene_id originale.\n"
  )
}


annotate_gene_ids <- function(gene_ids) {
  
  result <- data.frame(
    gene_id = gene_ids,
    Ensembl_ID = sub("\\..*$", "", gene_ids),
    gene_symbol = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (!annotation_available) {
    return(result)
  }
  
  unique_ids <- unique(result$Ensembl_ID)
  
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "list"
  )
  
  symbols <- vapply(
    result$Ensembl_ID,
    function(id) {
      
      value <- mapped[[id]]
      
      if (is.null(value) || length(value) == 0) {
        return(NA_character_)
      }
      
      value <- unique(
        value[
          !is.na(value) &
            value != ""
        ]
      )
      
      if (length(value) == 0) {
        return(NA_character_)
      }
      
      paste(
        value,
        collapse = ";"
      )
    },
    character(1)
  )
  
  result$gene_symbol <- symbols
  
  result
}


# ------------------------------------------------------------
# 5. ANALISI DEI CONTRASTI
# ------------------------------------------------------------

all_gene_results <- list()
all_category_summaries <- list()
all_method_summaries <- list()
all_correlation_summaries <- list()
all_profile_summaries <- list()

for (contrast in contrasts) {
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("CONTRASTO:", contrast, "\n")
  cat("------------------------------------------------------------\n")
  
  de_results <- lapply(
    methods,
    function(method) {
      read_de_result(
        method,
        contrast
      )
    }
  )
  
  names(de_results) <- methods
  
  reference_ids <- de_results[["TMM"]]$gene_id
  
  for (method in methods) {
    
    current_ids <- de_results[[method]]$gene_id
    
    if (!setequal(reference_ids, current_ids)) {
      stop(
        "Gene ID non allineati tra TMM e ",
        method,
        " nel contrasto ",
        contrast
      )
    }
    
    cat(
      method,
      ":",
      length(current_ids),
      "geni\n"
    )
  }
  
  gene_table <- data.frame(
    gene_id = reference_ids,
    stringsAsFactors = FALSE
  )
  
  for (method in methods) {
    
    current <- de_results[[method]]
    
    index <- match(
      reference_ids,
      current$gene_id
    )
    
    gene_table[[paste0("logFC_", method)]] <-
      current$logFC[index]
    
    gene_table[[paste0("FDR_", method)]] <-
      current$FDR[index]
  }
  
  
  # ----------------------------------------------------------
  # 5A. MATRICE DI SIGNIFICATIVITA'
  # ----------------------------------------------------------
  
  fdr_columns <- paste0(
    "FDR_",
    methods
  )
  
  logfc_columns <- paste0(
    "logFC_",
    methods
  )
  
  fdr_matrix <- as.matrix(
    gene_table[, fdr_columns, drop = FALSE]
  )
  
  logfc_matrix <- as.matrix(
    gene_table[, logfc_columns, drop = FALSE]
  )
  
  colnames(fdr_matrix) <- methods
  colnames(logfc_matrix) <- methods
  
  significance_matrix <- (
    !is.na(fdr_matrix) &
      fdr_matrix < fdr_threshold
  )
  
  colnames(significance_matrix) <- methods
  
  gene_table$significant_method_count <- rowSums(
    significance_matrix
  )
  
  gene_table$significant_methods <- apply(
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
  
  gene_table$significance_profile <- gene_table$significant_methods
  
  
  # ----------------------------------------------------------
  # 5B. CLASSIFICAZIONE PRINCIPALE
  # ----------------------------------------------------------
  
  tmm_sig <- significance_matrix[, "TMM"]
  
  alternative_sig <- significance_matrix[
    ,
    c("TMMwsp", "RLE", "upperquartile"),
    drop = FALSE
  ]
  
  any_alternative_sig <- rowSums(
    alternative_sig
  ) > 0
  
  all_alternative_sig <- rowSums(
    alternative_sig
  ) == ncol(alternative_sig)
  
  gene_table$discordance_class <- ifelse(
    tmm_sig & !any_alternative_sig,
    "TMM_ONLY",
    ifelse(
      !tmm_sig & any_alternative_sig,
      "ALTERNATIVE_ONLY",
      ifelse(
        tmm_sig & all_alternative_sig,
        "SIGNIFICANT_ALL_METHODS",
        ifelse(
          gene_table$significant_method_count == 0,
          "NOT_SIGNIFICANT_ANY_METHOD",
          "MIXED_METHOD_DISCORDANCE"
        )
      )
    )
  )
  
  
  # ----------------------------------------------------------
  # 5C. METRICHE FDR
  # ----------------------------------------------------------
  
  gene_table$minimum_FDR <- apply(
    fdr_matrix,
    1,
    safe_min
  )
  
  gene_table$maximum_FDR <- apply(
    fdr_matrix,
    1,
    safe_max
  )
  
  gene_table$FDR_range <- apply(
    fdr_matrix,
    1,
    safe_range
  )
  
  gene_table$distance_FDR_TMM_from_threshold <- abs(
    gene_table$FDR_TMM - fdr_threshold
  )
  
  gene_table$distance_minimum_FDR_from_threshold <- abs(
    gene_table$minimum_FDR - fdr_threshold
  )
  
  gene_table$near_FDR_threshold <- apply(
    fdr_matrix,
    1,
    function(x) {
      
      valid <- is.finite(x)
      
      any(
        valid &
          abs(x - fdr_threshold) <= near_fdr_window
      )
    }
  )
  
  gene_table$TMM_FDR_difference_max <- apply(
    abs(
      sweep(
        fdr_matrix[
          ,
          c("TMMwsp", "RLE", "upperquartile"),
          drop = FALSE
        ],
        1,
        fdr_matrix[, "TMM"],
        FUN = "-"
      )
    ),
    1,
    safe_max
  )
  
  
  # ----------------------------------------------------------
  # 5D. METRICHE LOGFC
  # ----------------------------------------------------------
  
  gene_table$minimum_logFC <- apply(
    logfc_matrix,
    1,
    safe_min
  )
  
  gene_table$maximum_logFC <- apply(
    logfc_matrix,
    1,
    safe_max
  )
  
  gene_table$logFC_range <- apply(
    logfc_matrix,
    1,
    safe_range
  )
  
  gene_table$max_abs_logFC_difference_vs_TMM <- apply(
    abs(
      sweep(
        logfc_matrix[
          ,
          c("TMMwsp", "RLE", "upperquartile"),
          drop = FALSE
        ],
        1,
        logfc_matrix[, "TMM"],
        FUN = "-"
      )
    ),
    1,
    safe_max
  )
  
  gene_table$direction_concordant_all <- apply(
    logfc_matrix,
    1,
    function(x) {
      length(unique(sign(x))) == 1
    }
  )
  
  gene_table$direction_concordant_vs_TMM <- apply(
    logfc_matrix[
      ,
      c("TMMwsp", "RLE", "upperquartile"),
      drop = FALSE
    ],
    1,
    function(x) {
      all(
        sign(x) == sign(
          logfc_matrix[, "TMM"][
            match(
              rownames(logfc_matrix),
              rownames(logfc_matrix)
            )
          ]
        )
      )
    }
  )
  
  # Calcolo esplicito e vettoriale della concordanza
  # rispetto alla direzione TMM.
  gene_table$direction_concordant_vs_TMM <- apply(
    logfc_matrix[
      ,
      c("TMMwsp", "RLE", "upperquartile"),
      drop = FALSE
    ],
    1,
    function(x) {
      TRUE
    }
  )
  
  gene_table$direction_concordant_vs_TMM <- (
    sign(logfc_matrix[, "TMM"]) ==
      sign(logfc_matrix[, "TMMwsp"]) &
      sign(logfc_matrix[, "TMM"]) ==
      sign(logfc_matrix[, "RLE"]) &
      sign(logfc_matrix[, "TMM"]) ==
      sign(logfc_matrix[, "upperquartile"])
  )
  
  
  # ----------------------------------------------------------
  # 5E. SIGNIFICATIVITA' ESCLUSIVA PER METODO
  # ----------------------------------------------------------
  
  for (method in methods) {
    
    other_methods <- setdiff(
      methods,
      method
    )
    
    other_significance <- significance_matrix[
      ,
      other_methods,
      drop = FALSE
    ]
    
    gene_table[[paste0(
      "significant_only_",
      method
    )]] <- (
      significance_matrix[, method] &
        rowSums(other_significance) == 0
    )
  }
  
  
  # ----------------------------------------------------------
  # 5F. ANNOTAZIONE GENI
  # ----------------------------------------------------------
  
  annotation <- annotate_gene_ids(
    gene_table$gene_id
  )
  
  gene_table <- merge(
    annotation,
    gene_table,
    by = "gene_id",
    all.y = TRUE,
    sort = FALSE
  )
  
  # Ripristino ordine originale dei gene ID.
  gene_table <- gene_table[
    match(reference_ids, gene_table$gene_id),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # 5G. ORDINAMENTO E OUTPUT COMPLETO
  # ----------------------------------------------------------
  
  gene_table <- gene_table[
    order(
      gene_table$discordance_class,
      gene_table$minimum_FDR,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  write.csv(
    gene_table,
    file.path(
      output_dir,
      paste0(
        "discordant_gene_audit_",
        contrast,
        ".csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  
  all_gene_results[[contrast]] <- gene_table
  
  
  # ----------------------------------------------------------
  # 5H. TABELLE PER CATEGORIA
  # ----------------------------------------------------------
  
  categories <- c(
    "TMM_ONLY",
    "ALTERNATIVE_ONLY",
    "SIGNIFICANT_ALL_METHODS",
    "MIXED_METHOD_DISCORDANCE",
    "NOT_SIGNIFICANT_ANY_METHOD"
  )
  
  category_summary <- data.frame(
    contrast = contrast,
    discordance_class = categories,
    gene_count = vapply(
      categories,
      function(category) {
        sum(
          gene_table$discordance_class == category,
          na.rm = TRUE
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  
  category_summary$percentage <- (
    category_summary$gene_count /
      nrow(gene_table)
  ) * 100
  
  all_category_summaries[[contrast]] <-
    category_summary
  
  for (category in categories) {
    
    subset_data <- gene_table[
      gene_table$discordance_class == category,
      ,
      drop = FALSE
    ]
    
    write.csv(
      subset_data,
      file.path(
        output_dir,
        paste0(
          "genes_",
          category,
          "_",
          contrast,
          ".csv"
        )
      ),
      row.names = FALSE,
      na = ""
    )
  }
  
  
  # ----------------------------------------------------------
  # 5I. GENI SIGNIFICATIVI ESCLUSIVI PER METODO
  # ----------------------------------------------------------
  
  method_specific_summary <- data.frame(
    contrast = contrast,
    method = methods,
    exclusive_significant_genes = vapply(
      methods,
      function(method) {
        sum(
          gene_table[[paste0(
            "significant_only_",
            method
          )]],
          na.rm = TRUE
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  
  all_method_summaries[[contrast]] <-
    method_specific_summary
  
  for (method in methods) {
    
    exclusive_genes <- gene_table[
      gene_table[[paste0(
        "significant_only_",
        method
      )]],
      ,
      drop = FALSE
    ]
    
    write.csv(
      exclusive_genes,
      file.path(
        output_dir,
        paste0(
          "genes_exclusive_",
          method,
          "_",
          contrast,
          ".csv"
        )
      ),
      row.names = FALSE,
      na = ""
    )
  }
  
  
  # ----------------------------------------------------------
  # 5J. GENI VICINI ALLA SOGLIA FDR
  # ----------------------------------------------------------
  
  near_threshold <- gene_table[
    gene_table$near_FDR_threshold,
    ,
    drop = FALSE
  ]
  
  near_threshold <- near_threshold[
    order(
      near_threshold$distance_FDR_TMM_from_threshold,
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
        contrast,
        ".csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  
  
  # ----------------------------------------------------------
  # 5K. GENI CON MAGGIORE VARIABILITA' LOGFC
  # ----------------------------------------------------------
  
  variable_logfc <- gene_table[
    gene_table$logFC_range >= 0.1,
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
        contrast,
        ".csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  
  
  # ----------------------------------------------------------
  # 5L. CONFRONTO TMM VS METODI ALTERNATIVI
  # ----------------------------------------------------------
  
  correlation_summary <- data.frame(
    contrast = contrast,
    reference = "TMM",
    alternative = setdiff(methods, "TMM"),
    Spearman_logFC = vapply(
      setdiff(methods, "TMM"),
      function(method) {
        safe_cor(
          logfc_matrix[, "TMM"],
          logfc_matrix[, method]
        )
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  all_correlation_summaries[[contrast]] <-
    correlation_summary
  
  
  # ----------------------------------------------------------
  # 5M. PROFILI DI SIGNIFICATIVITA'
  # ----------------------------------------------------------
  
  profile_summary <- as.data.frame(
    table(
      gene_table$significance_profile
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(profile_summary) <- c(
    "significance_profile",
    "gene_count"
  )
  
  profile_summary$contrast <- contrast
  
  profile_summary <- profile_summary[
    order(
      profile_summary$gene_count,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  all_profile_summaries[[contrast]] <-
    profile_summary
  
  write.csv(
    profile_summary,
    file.path(
      output_dir,
      paste0(
        "significance_profiles_",
        contrast,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 5N. REPORT A SCHERMO
  # ----------------------------------------------------------
  
  cat("\nClassificazione geni:\n")
  print(
    category_summary,
    row.names = FALSE
  )
  
  cat("\nGeni significativi esclusivi per metodo:\n")
  print(
    method_specific_summary,
    row.names = FALSE
  )
  
  cat("\nCorrelazioni logFC rispetto a TMM:\n")
  print(
    correlation_summary,
    row.names = FALSE
  )
  
  cat(
    "\nGeni vicini alla soglia FDR:",
    nrow(near_threshold),
    "\n"
  )
  
  cat(
    "Geni con range logFC >= 0.1:",
    nrow(variable_logfc),
    "\n"
  )
  
  cat(
    "Concordanza direzione logFC:",
    sum(
      gene_table$direction_concordant_all,
      na.rm = TRUE
    ),
    "/",
    nrow(gene_table),
    "\n"
  )
}


# ------------------------------------------------------------
# 6. RIEPILOGHI GLOBALI
# ------------------------------------------------------------

category_summary_all <- do.call(
  rbind,
  all_category_summaries
)

method_summary_all <- do.call(
  rbind,
  all_method_summaries
)

correlation_summary_all <- do.call(
  rbind,
  all_correlation_summaries
)

profile_summary_all <- do.call(
  rbind,
  all_profile_summaries
)

write.csv(
  category_summary_all,
  file.path(
    output_dir,
    "global_discordance_category_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  method_summary_all,
  file.path(
    output_dir,
    "global_method_exclusive_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  correlation_summary_all,
  file.path(
    output_dir,
    "global_logFC_correlation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  profile_summary_all,
  file.path(
    output_dir,
    "global_significance_profile_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. REPORT METODOLOGICO AUTOMATICO
# ------------------------------------------------------------

report_path <- file.path(
  output_dir,
  "NORMALIZATION_DISCORDANT_GENES_REPORT.txt"
)

sink(report_path)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("NORMALIZATION DISCORDANT GENES - FINAL REPORT\n")
cat("============================================================\n\n")

cat("Soglia FDR:", fdr_threshold, "\n")
cat("Intervallo near-threshold: +/-", near_fdr_window, "\n")
cat("Metodi:", paste(methods, collapse = ", "), "\n\n")

cat("CLASSIFICAZIONE GLOBALE\n")
cat("------------------------------------------------------------\n")
print(
  category_summary_all,
  row.names = FALSE
)

cat("\nGENI SIGNIFICATIVI ESCLUSIVI PER METODO\n")
cat("------------------------------------------------------------\n")
print(
  method_summary_all,
  row.names = FALSE
)

cat("\nCORRELAZIONI LOGFC\n")
cat("------------------------------------------------------------\n")
print(
  correlation_summary_all,
  row.names = FALSE
)

cat("\nPROFILI DI SIGNIFICATIVITA'\n")
cat("------------------------------------------------------------\n")
print(
  profile_summary_all,
  row.names = FALSE
)

cat("\n")
cat("INTERPRETAZIONE\n")
cat("------------------------------------------------------------\n")

cat(
  "TMM_ONLY: geni significativi con TMM ma non con gli altri metodi.\n"
)

cat(
  "ALTERNATIVE_ONLY: geni non significativi con TMM ma significativi con almeno un metodo alternativo.\n"
)

cat(
  "SIGNIFICANT_ALL_METHODS: geni significativi con tutti i metodi.\n"
)

cat(
  "MIXED_METHOD_DISCORDANCE: geni significativi con alcuni, ma non tutti, i metodi.\n"
)

cat(
  "NOT_SIGNIFICANT_ANY_METHOD: geni non significativi con nessun metodo.\n"
)

cat("\nLIMITI INTERPRETATIVI\n")
cat("------------------------------------------------------------\n")

cat(
  "La concordanza tra metodi di normalizzazione non dimostra l'assenza di bias.\n"
)

cat(
  "La significativita' dipendente dal metodo puo' riflettere differenze di FDR e precisione statistica.\n"
)

cat(
  "Questa analisi non sostituisce la disponibilita' di effective length RSEM originali.\n"
)

cat(
  "Le categorie discordanti non devono essere interpretate automaticamente come evidenza di effetti biologici distinti.\n"
)

sink()


# ------------------------------------------------------------
# 8. SALVATAGGIO OGGETTO COMPLETO
# ------------------------------------------------------------

saveRDS(
  list(
    gene_results = all_gene_results,
    category_summary = category_summary_all,
    method_exclusive_summary = method_summary_all,
    correlation_summary = correlation_summary_all,
    significance_profiles = profile_summary_all,
    configuration = list(
      methods = methods,
      contrasts = contrasts,
      fdr_threshold = fdr_threshold,
      near_fdr_window = near_fdr_window
    ),
    annotation_available = annotation_available,
    session_info = sessionInfo()
  ),
  file = file.path(
    output_dir,
    "normalization_discordant_genes_complete.rds"
  )
)


# ------------------------------------------------------------
# 9. SESSION INFO
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 10. INVENTARIO OUTPUT
# ------------------------------------------------------------

generated_files <- list.files(
  output_dir,
  full.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("NORMALIZATION DISCORDANT GENES - FINAL REPORT\n")
cat("============================================================\n\n")

cat("Riepilogo classificazione:\n")
print(
  category_summary_all,
  row.names = FALSE
)

cat("\nGeni significativi esclusivi per metodo:\n")
print(
  method_summary_all,
  row.names = FALSE
)

cat("\nCorrelazioni logFC:\n")
print(
  correlation_summary_all,
  row.names = FALSE
)

cat("\nDirectory output:\n")
cat(normalizePath(output_dir), "\n")

cat("\nNumero file generati:", length(generated_files), "\n")

cat("\nFile generati:\n")
print(generated_files)

cat("\n")
cat("ANALISI COMPLETATA.\n")
