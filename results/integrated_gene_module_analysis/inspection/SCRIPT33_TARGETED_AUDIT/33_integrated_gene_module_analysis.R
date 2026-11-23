# ============================================================
# TRANSCRIPTOME-X
# Script: 33_integrated_gene_module_analysis.R
#
# INTEGRATED GENE MODULE ANALYSIS
#
# Analisi descrittiva di quattro moduli genici candidati:
# 1. ECM remodeling
# 2. TGF-beta/BMP signaling
# 3. Cell adhesion/cytoskeleton
# 4. Growth factor/receptor signaling
#
# Integra i due contrasti DE TMM, confronta logFC/FDR/direzione,
# produce audit, tabelle, figure e report.
#
# I moduli sono liste candidate definite a priori:
# NON costituiscono un test formale di pathway enrichment e
# NON dimostrano attivazione di pathway o causalita'.
#
# Requisiti: R >= 4.0; solo pacchetti base R.
# ============================================================


# ============================================================
# 0. CONFIGURAZIONE
# ============================================================

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

PROJECT_ROOT <- "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"

INPUT_DE_DIR <- file.path(
  PROJECT_ROOT, "results", "normalization_sensitivity"
)

INPUT_SIGNATURES <- file.path(
  PROJECT_ROOT,
  "results",
  "functional_enrichment_interpretation",
  "tables",
  "significant_term_gene_signatures.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT, "results", "integrated_gene_module_analysis"
)

OUTPUT_TABLES <- file.path(OUTPUT_DIR, "tables")
OUTPUT_FIGURES <- file.path(OUTPUT_DIR, "figures")
OUTPUT_AUDIT <- file.path(OUTPUT_DIR, "audit")
OUTPUT_REPORTS <- file.path(OUTPUT_DIR, "reports")

FDR_CUTOFF <- 0.05

CONTRAST_1 <- "siNRP1_3_vs_siNT"
CONTRAST_2 <- "siNRP1_5_vs_siNT"

# Nomi effettivi dei file DE condivisi in precedenza.
FILE_DE_1 <- "DE_TMM_siNRP1_3_vs_siNT.csv"
FILE_DE_2 <- "DE_TMM_siNRP1_5_vs_siNT.csv"

SIGNATURE_FILE <- "significant_term_gene_signatures.csv"


# ============================================================
# 1. MODULI GENICI CANDIDATI
# ============================================================

GENE_MODULES <- list(

  ECM_REMODELING = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
    "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "FN1", "FBLN1", "FBLN2", "FBLN5", "LUM", "DCN",
    "VCAN", "BGN", "THBS1", "THBS2", "THBS3",
    "MMP1", "MMP2", "MMP3", "MMP9", "MMP14",
    "TIMP1", "TIMP2", "TIMP3", "TIMP4",
    "ADAMTS1", "ADAMTS4", "ADAMTS5", "ADAMTS6",
    "PLOD1", "PLOD2", "LOX", "LOXL1", "LOXL2",
    "HTRA1", "MATN2", "MATN3", "FRAS1",
    "CCN2", "CCN3", "RELN", "VTN", "TNC",
    "POSTN", "SPARC", "SERPINE1", "SERPINH1"
  ),

  TGFB_BMP_SIGNALING = c(
    "TGFB1", "TGFB2", "TGFB3",
    "TGFBR1", "TGFBR2", "TGFBR3",
    "SMAD2", "SMAD3", "SMAD4", "SMAD6", "SMAD7",
    "BMP2", "BMP4", "BMP5", "BMP6", "BMP7",
    "BMPR1A", "BMPR1B", "BMPR2",
    "NOG", "GREM1", "GREM2",
    "ID1", "ID2", "ID3",
    "CCN2", "CCN3", "THBS1",
    "ENG", "ACVR1", "ACVR2A", "ACVR2B"
  ),

  CELL_ADHESION_CYTOSKELETON = c(
    "CAV1", "CAV2", "ITGA1", "ITGA2", "ITGA3",
    "ITGA5", "ITGA6", "ITGAV",
    "ITGB1", "ITGB3", "ITGB4", "ITGB5",
    "PXN", "VCL", "TLN1", "TLN2", "FERMT2",
    "ACTN1", "ACTN4", "FLNA", "FLNB",
    "VIM", "TAGLN", "TPM1", "TPM2", "TPM3",
    "MYH9", "MYH10", "MYL9", "CNN1",
    "RHOA", "RHOC", "RAC1", "CDC42",
    "ROCK1", "ROCK2", "PAK1", "PAK2",
    "GPC1", "GPC4", "SDC1", "SDC2", "SDC4",
    "CTNNB1", "CDH1", "CDH2",
    "PDLIM1", "PDLIM5", "FSCN1"
  ),

  GROWTH_FACTOR_RECEPTOR_SIGNALING = c(
    "NRP1", "VEGFA", "VEGFB", "VEGFC",
    "FLT1", "KDR", "FLT4",
    "FGF1", "FGF2", "FGF7", "FGF10",
    "FGFR1", "FGFR2", "FGFR3",
    "EGF", "EGFR", "ERBB2", "ERBB3",
    "IGF1", "IGF2", "IGF1R",
    "IGFBP3", "IGFBP5",
    "MET", "HGF",
    "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
    "PIK3CA", "PIK3CB", "AKT1", "AKT2",
    "GRB2", "SOS1", "KRAS", "BRAF", "MAPK1",
    "MAPK3", "IL6", "IL6R", "IL6ST",
    "LYN", "SRC", "CCN2"
  )
)


# ============================================================
# 2. FUNZIONI GENERALI
# ============================================================

log_message <- function(...) {
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "|",
    paste(..., collapse = " "),
    "\n"
  )
}

stop_with_message <- function(...) {
  stop(paste(..., collapse = " "), call. = FALSE)
}

write_csv_safe <- function(data, path) {
  if (is.null(data)) {
    data <- data.frame()
  }

  write.csv(
    data,
    file = path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

normalize_column_names <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]", "", x))
}

find_column <- function(data, candidates, required = TRUE) {
  original_names <- names(data)
  normalized_names <- normalize_column_names(original_names)
  normalized_candidates <- normalize_column_names(candidates)

  index <- match(normalized_candidates, normalized_names)
  index <- index[!is.na(index)]

  if (length(index) > 0) {
    return(original_names[index[1]])
  }

  if (required) {
    stop_with_message(
      "Colonna obbligatoria non trovata.",
      "Candidate:", paste(candidates, collapse = ", "),
      "Colonne disponibili:", paste(original_names, collapse = ", ")
    )
  }

  NULL
}

# Ricerca letterale del nome file: evita problemi di escaping regex.
find_input_file <- function(root, filename, preferred_dir = NULL) {
  if (!dir.exists(root)) {
    stop_with_message("Directory inesistente:", root)
  }

  if (!is.null(preferred_dir) && dir.exists(preferred_dir)) {
    preferred_path <- file.path(preferred_dir, filename)

    if (file.exists(preferred_path)) {
      return(normalizePath(
        preferred_path, winslash = "/", mustWork = TRUE
      ))
    }
  }

  matches <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )

  matches <- matches[
    tolower(basename(matches)) == tolower(filename) &
      file.info(matches)$isdir %in% FALSE
  ]

  if (length(matches) == 0) {
    stop_with_message(
      "File non trovato:", filename,
      "nella directory:", root
    )
  }

  if (length(matches) > 1) {
    log_message(
      "ATTENZIONE: trovate", length(matches),
      "copie di", filename,
      "- seleziono quella nella directory preferita, se presente."
    )

    if (!is.null(preferred_dir)) {
      preferred_matches <- matches[
        startsWith(
          tolower(normalizePath(matches, winslash = "/", mustWork = TRUE)),
          tolower(normalizePath(preferred_dir, winslash = "/", mustWork = TRUE))
        )
      ]

      if (length(preferred_matches) > 0) {
        matches <- preferred_matches
      }
    }
  }

  normalizePath(matches[1], winslash = "/", mustWork = TRUE)
}

read_csv_checked <- function(path) {
  if (!file.exists(path)) {
    stop_with_message("File inesistente:", path)
  }

  data <- tryCatch(
    read.csv(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "NaN")
    ),
    error = function(e) {
      stop_with_message(
        "Errore durante la lettura di", basename(path),
        ":", conditionMessage(e)
      )
    }
  )

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop_with_message("File CSV vuoto:", path)
  }

  data
}

as_numeric_checked <- function(x, column_name) {
  original <- trimws(as.character(x))
  converted <- suppressWarnings(as.numeric(original))

  invalid <- !is.na(original) & nzchar(original) & is.na(converted)

  if (any(invalid)) {
    stop_with_message(
      "Valori non numerici nella colonna", column_name,
      "- valori problematici:", sum(invalid)
    )
  }

  converted
}

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_
  x
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_
  x
}

safe_sum <- function(x) {
  sum(x, na.rm = TRUE)
}


# ============================================================
# 3. PREPARAZIONE DIRECTORY E INPUT
# ============================================================

if (!dir.exists(PROJECT_ROOT)) {
  stop_with_message("PROJECT_ROOT non trovato:", PROJECT_ROOT)
}

dir.create(OUTPUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_AUDIT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_REPORTS, recursive = TRUE, showWarnings = FALSE)

log_message("TRANSCRIPTOME-X - SCRIPT 33")
log_message("Project root:", PROJECT_ROOT)
log_message("Output:", OUTPUT_DIR)

de_file_1 <- find_input_file(PROJECT_ROOT, FILE_DE_1, INPUT_DE_DIR)
de_file_2 <- find_input_file(PROJECT_ROOT, FILE_DE_2, INPUT_DE_DIR)
signature_path <- find_input_file(
  PROJECT_ROOT, SIGNATURE_FILE, dirname(INPUT_SIGNATURES)
)

log_message("DE file 1:", de_file_1)
log_message("DE file 2:", de_file_2)
log_message("Signatures:", signature_path)

de_raw_1 <- read_csv_checked(de_file_1)
de_raw_2 <- read_csv_checked(de_file_2)
signatures_raw <- read_csv_checked(signature_path)


# ============================================================
# 4. STANDARDIZZAZIONE RISULTATI DE
# ============================================================

standardize_de <- function(data, contrast_name) {
  symbol_col <- find_column(
    data,
    c(
      "gene_symbol", "symbol", "gene_name",
      "external_gene_name", "genesymbol", "Gene",
      "hgnc_symbol", "external_gene",
      "Ensembl.114.Gene.Name"
    ),
    required = TRUE
  )

  logfc_col <- find_column(
    data,
    c("logFC", "log2FoldChange", "log_fold_change"),
    required = TRUE
  )

  fdr_col <- find_column(
    data,
    c("FDR", "padj", "adjusted_p_value", "p_adj", "adj.P.Val"),
    required = TRUE
  )

  id_col <- find_column(
    data,
    c(
      "gene_id", "ensembl_gene_id", "ensembl_id",
      "GeneID", "gene", "ID", "Geneid",
      "Ensembl.114.Gene.ID"
    ),
    required = FALSE
  )

  pvalue_col <- find_column(
    data, c("PValue", "pvalue", "p_value", "P.Value"),
    required = FALSE
  )

  logcpm_col <- find_column(
    data, c("logCPM", "log_cpm"),
    required = FALSE
  )

  result <- data.frame(
    contrast = contrast_name,
    gene_symbol = clean_symbol(data[[symbol_col]]),
    gene_id = if (!is.null(id_col)) {
      clean_gene_id(data[[id_col]])
    } else {
      rep(NA_character_, nrow(data))
    },
    logFC = as_numeric_checked(data[[logfc_col]], logfc_col),
    FDR = as_numeric_checked(data[[fdr_col]], fdr_col),
    stringsAsFactors = FALSE
  )

  result$PValue <- if (!is.null(pvalue_col)) {
    as_numeric_checked(data[[pvalue_col]], pvalue_col)
  } else {
    rep(NA_real_, nrow(data))
  }

  result$logCPM <- if (!is.null(logcpm_col)) {
    as_numeric_checked(data[[logcpm_col]], logcpm_col)
  } else {
    rep(NA_real_, nrow(data))
  }

  if (anyNA(result$logFC)) {
    stop_with_message("logFC mancanti nel contrasto:", contrast_name)
  }

  if (anyNA(result$FDR)) {
    stop_with_message("FDR mancanti nel contrasto:", contrast_name)
  }

  if (any(result$FDR < 0 | result$FDR > 1)) {
    stop_with_message("FDR fuori da [0,1] nel contrasto:", contrast_name)
  }

  result$significant <- result$FDR < FDR_CUTOFF
  result$direction <- ifelse(
    result$logFC > 0, "UP",
    ifelse(result$logFC < 0, "DOWN", "ZERO")
  )

  result
}

de_1_all <- standardize_de(de_raw_1, CONTRAST_1)
de_2_all <- standardize_de(de_raw_2, CONTRAST_2)


# ============================================================
# 5. AUDIT DUPLICATI E RISOLUZIONE SIMBOLI
# ============================================================

audit_duplicate_symbols <- function(data) {
  valid <- data[!is.na(data$gene_symbol), , drop = FALSE]

  if (nrow(valid) == 0) {
    return(data.frame(
      gene_symbol = character(),
      number_of_records = integer()
    ))
  }

  counts <- as.data.frame(
    table(valid$gene_symbol),
    stringsAsFactors = FALSE
  )
  names(counts) <- c("gene_symbol", "number_of_records")
  counts <- counts[counts$number_of_records > 1, , drop = FALSE]
  counts$number_of_records <- as.integer(counts$number_of_records)
  rownames(counts) <- NULL
  counts
}

duplicate_audit_1 <- audit_duplicate_symbols(de_1_all)
duplicate_audit_2 <- audit_duplicate_symbols(de_2_all)

if (nrow(duplicate_audit_1) > 0) {
  duplicate_audit_1$contrast <- CONTRAST_1
}
if (nrow(duplicate_audit_2) > 0) {
  duplicate_audit_2$contrast <- CONTRAST_2
}

duplicate_audit <- rbind(duplicate_audit_1, duplicate_audit_2)

write_csv_safe(
  duplicate_audit,
  file.path(OUTPUT_AUDIT, "duplicate_gene_symbol_audit.csv")
)

# Per il confronto tra contrasti viene mantenuto un record per simbolo:
# priorita' al FDR piu' basso; a parita' di FDR, maggiore |logFC|.
# Tutti i duplicati originali restano documentati nell'audit.
select_best_symbol_record <- function(data) {
  data <- data[!is.na(data$gene_symbol), , drop = FALSE]

  data <- data[
    order(
      data$gene_symbol,
      data$FDR,
      -abs(data$logFC),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  data <- data[!duplicated(data$gene_symbol), , drop = FALSE]
  rownames(data) <- NULL
  data
}

de_1 <- select_best_symbol_record(de_1_all)
de_2 <- select_best_symbol_record(de_2_all)

log_message("DE records:", nrow(de_1_all), "and", nrow(de_2_all))
log_message(
  "Unique gene symbols:",
  length(unique(de_1$gene_symbol)), "and",
  length(unique(de_2$gene_symbol))
)


# ============================================================
# 6. AUDIT INPUT E CONFRONTO GENERALE
# ============================================================

input_audit <- data.frame(
  contrast = c(CONTRAST_1, CONTRAST_2),
  input_file = c(de_file_1, de_file_2),
  total_records = c(nrow(de_1_all), nrow(de_2_all)),
  records_with_symbol = c(
    sum(!is.na(de_1_all$gene_symbol)),
    sum(!is.na(de_2_all$gene_symbol))
  ),
  unique_symbols = c(
    length(unique(na.omit(de_1_all$gene_symbol))),
    length(unique(na.omit(de_2_all$gene_symbol)))
  ),
  significant_genes = c(
    sum(de_1_all$significant),
    sum(de_2_all$significant)
  ),
  up_genes = c(
    sum(de_1_all$significant & de_1_all$direction == "UP"),
    sum(de_2_all$significant & de_2_all$direction == "UP")
  ),
  down_genes = c(
    sum(de_1_all$significant & de_1_all$direction == "DOWN"),
    sum(de_2_all$significant & de_2_all$direction == "DOWN")
  )
)

write_csv_safe(
  input_audit,
  file.path(OUTPUT_AUDIT, "input_audit.csv")
)

de_comparison <- merge(
  de_1,
  de_2,
  by = "gene_symbol",
  suffixes = c("_siNRP1_3", "_siNRP1_5"),
  all = TRUE
)

de_comparison$direction_concordance <- ifelse(
  is.na(de_comparison$logFC_siNRP1_3) |
    is.na(de_comparison$logFC_siNRP1_5),
  "NOT_COMPARABLE",
  ifelse(
    sign(de_comparison$logFC_siNRP1_3) ==
      sign(de_comparison$logFC_siNRP1_5),
    "CONCORDANT",
    "DISCORDANT"
  )
)

sig1 <- !is.na(de_comparison$significant_siNRP1_3) &
  de_comparison$significant_siNRP1_3
sig2 <- !is.na(de_comparison$significant_siNRP1_5) &
  de_comparison$significant_siNRP1_5

de_comparison$significant_both <- sig1 & sig2
de_comparison$significant_either <- sig1 | sig2

de_comparison$joint_class <- ifelse(
  de_comparison$significant_both &
    de_comparison$direction_concordance == "CONCORDANT",
  "SIGNIFICANT_BOTH_CONCORDANT",
  ifelse(
    de_comparison$significant_both &
      de_comparison$direction_concordance == "DISCORDANT",
    "SIGNIFICANT_BOTH_DISCORDANT",
    ifelse(
      sig1,
      "SIGNIFICANT_ONLY_siNRP1_3",
      ifelse(sig2, "SIGNIFICANT_ONLY_siNRP1_5", "NOT_SIGNIFICANT_BOTH")
    )
  )
)

write_csv_safe(
  de_comparison,
  file.path(OUTPUT_TABLES, "global_DE_gene_comparison.csv")
)


# ============================================================
# 7. DEFINIZIONE E PRESENZA DEI MODULI
# ============================================================

module_gene_table <- do.call(
  rbind,
  lapply(names(GENE_MODULES), function(module_name) {
    data.frame(
      module = module_name,
      gene_symbol = unique(GENE_MODULES[[module_name]]),
      stringsAsFactors = FALSE
    )
  })
)

rownames(module_gene_table) <- NULL

write_csv_safe(
  module_gene_table,
  file.path(OUTPUT_TABLES, "module_gene_definitions.csv")
)

all_module_genes <- unique(module_gene_table$gene_symbol)

module_presence <- merge(
  module_gene_table,
  data.frame(
    gene_symbol = all_module_genes,
    in_siNRP1_3 = all_module_genes %in% de_1$gene_symbol,
    in_siNRP1_5 = all_module_genes %in% de_2$gene_symbol,
    stringsAsFactors = FALSE
  ),
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

write_csv_safe(
  module_presence,
  file.path(OUTPUT_AUDIT, "module_gene_presence_audit.csv")
)


# ============================================================
# 8. RISULTATI DE PER MODULO
# ============================================================

extract_module_results <- function(module_name, genes, de_data) {
  selected <- de_data[
    de_data$gene_symbol %in% genes,
    ,
    drop = FALSE
  ]

  if (nrow(selected) == 0) {
    return(data.frame(
      module = character(),
      contrast = character(),
      gene_symbol = character(),
      gene_id = character(),
      logFC = numeric(),
      FDR = numeric(),
      significant = logical(),
      direction = character(),
      PValue = numeric(),
      logCPM = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  selected$module <- module_name

  selected[, c(
    "module", "contrast", "gene_symbol", "gene_id",
    "logFC", "FDR", "significant", "direction",
    "PValue", "logCPM"
  ), drop = FALSE]
}

module_results_list <- list()

for (module_name in names(GENE_MODULES)) {
  genes <- GENE_MODULES[[module_name]]

  module_results_list[[paste0(module_name, "_siNRP1_3")]] <-
    extract_module_results(module_name, genes, de_1)

  module_results_list[[paste0(module_name, "_siNRP1_5")]] <-
    extract_module_results(module_name, genes, de_2)
}

module_results <- do.call(rbind, module_results_list)
rownames(module_results) <- NULL

write_csv_safe(
  module_results,
  file.path(OUTPUT_TABLES, "module_gene_results_long.csv")
)


# ============================================================
# 9. CONFRONTO INTEGRATO DEI MODULI
# ============================================================

module_comparison_list <- list()

for (module_name in names(GENE_MODULES)) {
  genes <- GENE_MODULES[[module_name]]

  m1 <- de_1[de_1$gene_symbol %in% genes, , drop = FALSE]
  m2 <- de_2[de_2$gene_symbol %in% genes, , drop = FALSE]

  comparison <- merge(
    m1,
    m2,
    by = "gene_symbol",
    suffixes = c("_siNRP1_3", "_siNRP1_5"),
    all = TRUE,
    sort = TRUE
  )

  comparison$module <- module_name

  comparison$direction_concordance <- ifelse(
    is.na(comparison$logFC_siNRP1_3) |
      is.na(comparison$logFC_siNRP1_5),
    "NOT_COMPARABLE",
    ifelse(
      sign(comparison$logFC_siNRP1_3) ==
        sign(comparison$logFC_siNRP1_5),
      "CONCORDANT",
      "DISCORDANT"
    )
  )

  sig1_module <- !is.na(comparison$significant_siNRP1_3) &
    comparison$significant_siNRP1_3
  sig2_module <- !is.na(comparison$significant_siNRP1_5) &
    comparison$significant_siNRP1_5

  comparison$significant_both <- sig1_module & sig2_module
  comparison$significant_either <- sig1_module | sig2_module

  comparison$joint_class <- ifelse(
    comparison$significant_both &
      comparison$direction_concordance == "CONCORDANT",
    "SIGNIFICANT_BOTH_CONCORDANT",
    ifelse(
      comparison$significant_both &
        comparison$direction_concordance == "DISCORDANT",
      "SIGNIFICANT_BOTH_DISCORDANT",
      ifelse(
        sig1_module,
        "SIGNIFICANT_ONLY_siNRP1_3",
        ifelse(
          sig2_module,
          "SIGNIFICANT_ONLY_siNRP1_5",
          "NOT_SIGNIFICANT_BOTH"
        )
      )
    )
  )

  module_comparison_list[[module_name]] <- comparison
}

module_comparison <- do.call(rbind, module_comparison_list)
rownames(module_comparison) <- NULL

write_csv_safe(
  module_comparison,
  file.path(OUTPUT_TABLES, "module_integrated_gene_comparison.csv")
)


# ============================================================
# 10. RIEPILOGO MODULI
# ============================================================

module_summary <- do.call(
  rbind,
  lapply(names(GENE_MODULES), function(module_name) {
    genes <- GENE_MODULES[[module_name]]

    m1 <- de_1[de_1$gene_symbol %in% genes, , drop = FALSE]
    m2 <- de_2[de_2$gene_symbol %in% genes, , drop = FALSE]

    common_symbols <- intersect(m1$gene_symbol, m2$gene_symbol)

    both <- module_comparison[
      module_comparison$module == module_name &
        module_comparison$significant_both,
      ,
      drop = FALSE
    ]

    data.frame(
      module = module_name,
      candidate_genes = length(unique(genes)),
      detected_siNRP1_3 = length(unique(m1$gene_symbol)),
      detected_siNRP1_5 = length(unique(m2$gene_symbol)),
      detected_both = length(common_symbols),
      significant_siNRP1_3 = sum(m1$significant),
      significant_siNRP1_5 = sum(m2$significant),
      significant_both = nrow(both),
      significant_both_concordant = sum(
        both$direction_concordance == "CONCORDANT",
        na.rm = TRUE
      ),
      significant_both_discordant = sum(
        both$direction_concordance == "DISCORDANT",
        na.rm = TRUE
      ),
      up_siNRP1_3 = sum(m1$significant & m1$direction == "UP"),
      down_siNRP1_3 = sum(m1$significant & m1$direction == "DOWN"),
      up_siNRP1_5 = sum(m2$significant & m2$direction == "UP"),
      down_siNRP1_5 = sum(m2$significant & m2$direction == "DOWN"),
      stringsAsFactors = FALSE
    )
  })
)

write_csv_safe(
  module_summary,
  file.path(OUTPUT_TABLES, "module_summary.csv")
)


# ============================================================
# 11. GENI PRIORITARI
# ============================================================

module_priority_genes <- module_comparison[
  module_comparison$significant_either,
  ,
  drop = FALSE
]

if (nrow(module_priority_genes) > 0) {
  module_priority_genes$minimum_FDR <- pmin(
    module_priority_genes$FDR_siNRP1_3,
    module_priority_genes$FDR_siNRP1_5,
    na.rm = TRUE
  )

  module_priority_genes$maximum_absolute_logFC <- pmax(
    abs(module_priority_genes$logFC_siNRP1_3),
    abs(module_priority_genes$logFC_siNRP1_5),
    na.rm = TRUE
  )

  module_priority_genes <- module_priority_genes[
    order(
      module_priority_genes$module,
      module_priority_genes$minimum_FDR,
      -module_priority_genes$maximum_absolute_logFC,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}

write_csv_safe(
  module_priority_genes,
  file.path(OUTPUT_TABLES, "module_priority_genes.csv")
)


# ============================================================
# 12. INTEGRAZIONE SIGNATURE DI ENRICHMENT
# ============================================================

signature_gene_col <- find_column(
  signatures_raw,
  c(
    "gene_symbol", "gene_symbols", "symbols",
    "genesymbol", "gene_list", "genes", "geneid",
    "gene_ids", "geneID"
  ),
  required = FALSE
)

signature_term_col <- find_column(
  signatures_raw,
  c("term_description", "Description", "term", "pathway"),
  required = FALSE
)

signature_contrast_col <- find_column(
  signatures_raw,
  c("contrast", "comparison"),
  required = FALSE
)

signature_group_col <- find_column(
  signatures_raw,
  c("group", "gene_group", "method_group"),
  required = FALSE
)

signature_ontology_col <- find_column(
  signatures_raw,
  c("ontology", "ont", "category"),
  required = FALSE
)

signature_audit <- data.frame(
  file = signature_path,
  rows = nrow(signatures_raw),
  columns = ncol(signatures_raw),
  gene_column_detected = if (is.null(signature_gene_col)) {
    NA_character_
  } else {
    signature_gene_col
  },
  term_column_detected = if (is.null(signature_term_col)) {
    NA_character_
  } else {
    signature_term_col
  },
  contrast_column_detected = if (is.null(signature_contrast_col)) {
    NA_character_
  } else {
    signature_contrast_col
  },
  group_column_detected = if (is.null(signature_group_col)) {
    NA_character_
  } else {
    signature_group_col
  },
  ontology_column_detected = if (is.null(signature_ontology_col)) {
    NA_character_
  } else {
    signature_ontology_col
  },
  stringsAsFactors = FALSE
)

write_csv_safe(
  signature_audit,
  file.path(OUTPUT_AUDIT, "signature_input_audit.csv")
)

signature_module_overlap <- data.frame()
signature_module_summary <- data.frame()

if (is.null(signature_gene_col)) {
  log_message(
    "Colonna genica non riconosciuta nel file signatures.",
    "Analisi DE dei moduli completata comunque."
  )
} else {
  split_gene_tokens <- function(x) {
    if (is.na(x) || !nzchar(trimws(as.character(x)))) {
      return(character())
    }

    tokens <- unlist(
      strsplit(as.character(x), "[/;,|[:space:]]+"),
      use.names = FALSE
    )

    tokens <- clean_symbol(tokens)

    # Si conservano solo token alfabetici compatibili con simboli.
    # Gli identificativi numerici Entrez non vengono convertiti.
    tokens <- tokens[
      !is.na(tokens) &
        grepl("[A-Z]", tokens) &
        !grepl("^ENS[A-Z]*G[0-9]+", tokens)
    ]

    unique(tokens)
  }

  signature_rows <- vector("list", nrow(signatures_raw))

  for (i in seq_len(nrow(signatures_raw))) {
    genes <- split_gene_tokens(signatures_raw[[signature_gene_col]][i])

    if (length(genes) == 0) {
      next
    }

    signature_rows[[i]] <- data.frame(
      signature_row = i,
      contrast = if (is.null(signature_contrast_col)) {
        NA_character_
      } else {
        as.character(signatures_raw[[signature_contrast_col]][i])
      },
      group = if (is.null(signature_group_col)) {
        NA_character_
      } else {
        as.character(signatures_raw[[signature_group_col]][i])
      },
      ontology = if (is.null(signature_ontology_col)) {
        NA_character_
      } else {
        as.character(signatures_raw[[signature_ontology_col]][i])
      },
      term = if (is.null(signature_term_col)) {
        NA_character_
      } else {
        as.character(signatures_raw[[signature_term_col]][i])
      },
      gene_symbol = genes,
      stringsAsFactors = FALSE
    )
  }

  signature_rows <- Filter(Negate(is.null), signature_rows)

  if (length(signature_rows) == 0) {
    log_message(
      "Nessuna firma genica interpretabile come simboli genici."
    )
  } else {
    signature_gene_long <- do.call(rbind, signature_rows)
    rownames(signature_gene_long) <- NULL

    signature_module_overlap <- merge(
      signature_gene_long,
      module_gene_table,
      by = "gene_symbol",
      all = FALSE,
      sort = FALSE
    )

    signature_module_overlap <- unique(signature_module_overlap)

    write_csv_safe(
      signature_module_overlap,
      file.path(OUTPUT_TABLES, "signature_module_gene_overlap.csv")
    )

    if (nrow(signature_module_overlap) > 0) {
      signature_module_summary <- aggregate(
        gene_symbol ~ contrast + group + ontology + term + module,
        data = signature_module_overlap,
        FUN = function(x) length(unique(x))
      )

      names(signature_module_summary)[
        names(signature_module_summary) == "gene_symbol"
      ] <- "module_genes_in_signature"
    }

    write_csv_safe(
      signature_module_summary,
      file.path(OUTPUT_TABLES, "signature_module_summary.csv")
    )
  }
}


# ============================================================
# 13. FIGURE
# ============================================================

# 13.1 Conteggio geni significativi per modulo.
png(
  filename = file.path(
    OUTPUT_FIGURES, "module_significant_gene_counts.png"
  ),
  width = 1600,
  height = 1000,
  res = 150
)

par(mar = c(10, 5, 4, 2))

count_matrix <- rbind(
  siNRP1_3 = module_summary$significant_siNRP1_3,
  siNRP1_5 = module_summary$significant_siNRP1_5
)

barplot(
  count_matrix,
  beside = TRUE,
  names.arg = module_summary$module,
  las = 2,
  ylab = "Numero di geni significativi",
  main = "Significant genes per candidate module",
  legend.text = rownames(count_matrix),
  args.legend = list(x = "topright", bty = "n")
)

dev.off()


# 13.2 Scatter plot logFC per modulo.
for (module_name in names(GENE_MODULES)) {
  plot_data <- module_comparison[
    module_comparison$module == module_name &
      !is.na(module_comparison$logFC_siNRP1_3) &
      !is.na(module_comparison$logFC_siNRP1_5),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    next
  }

  png(
    filename = file.path(
      OUTPUT_FIGURES,
      paste0(module_name, "_logFC_scatter.png")
    ),
    width = 1200,
    height = 1000,
    res = 150
  )

  plot(
    plot_data$logFC_siNRP1_3,
    plot_data$logFC_siNRP1_5,
    pch = 19,
    cex = 0.8,
    xlab = "logFC siNRP1.3 vs siNT",
    ylab = "logFC siNRP1.5 vs siNT",
    main = module_name
  )

  abline(h = 0, v = 0, lty = 2)
  abline(a = 0, b = 1, lty = 3)

  concordant <- plot_data$direction_concordance == "CONCORDANT"

  points(
    plot_data$logFC_siNRP1_3[concordant],
    plot_data$logFC_siNRP1_5[concordant],
    pch = 19
  )

  points(
    plot_data$logFC_siNRP1_3[!concordant],
    plot_data$logFC_siNRP1_5[!concordant],
    pch = 1
  )

  legend(
    "topleft",
    legend = c("Direzione concordante", "Direzione discordante"),
    pch = c(19, 1),
    bty = "n"
  )

  dev.off()
}


# 13.3 Heatmap logFC.
heatmap_data <- module_comparison[
  !is.na(module_comparison$logFC_siNRP1_3) &
    !is.na(module_comparison$logFC_siNRP1_5),
  ,
  drop = FALSE
]

if (nrow(heatmap_data) > 0) {
  heatmap_matrix <- as.matrix(
    heatmap_data[, c("logFC_siNRP1_3", "logFC_siNRP1_5")]
  )

  storage.mode(heatmap_matrix) <- "numeric"

  rownames(heatmap_matrix) <- make.unique(
    paste(heatmap_data$module, heatmap_data$gene_symbol, sep = " | ")
  )

  png(
    filename = file.path(
      OUTPUT_FIGURES, "integrated_module_logFC_heatmap.png"
    ),
    width = 1100,
    height = max(1000, min(5000, nrow(heatmap_matrix) * 18)),
    res = 150
  )

  heatmap(
    heatmap_matrix,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(12, 18),
    main = "Integrated candidate-module gene logFC"
  )

  dev.off()
}


# ============================================================
# 14. AUDIT FINALE
# ============================================================

final_audit <- data.frame(
  metric = c(
    "DE records siNRP1.3",
    "DE records siNRP1.5",
    "Unique symbols siNRP1.3",
    "Unique symbols siNRP1.5",
    "Common symbols between contrasts",
    "Candidate module genes (unique)",
    "Module comparison records",
    "Module genes significant in both",
    "Module genes concordant in both",
    "Module genes discordant in both"
  ),
  value = c(
    nrow(de_1_all),
    nrow(de_2_all),
    length(unique(de_1$gene_symbol)),
    length(unique(de_2$gene_symbol)),
    nrow(de_comparison),
    length(all_module_genes),
    nrow(module_comparison),
    sum(module_comparison$significant_both, na.rm = TRUE),
    sum(
      module_comparison$significant_both &
        module_comparison$direction_concordance == "CONCORDANT",
      na.rm = TRUE
    ),
    sum(
      module_comparison$significant_both &
        module_comparison$direction_concordance == "DISCORDANT",
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  final_audit,
  file.path(OUTPUT_AUDIT, "final_analysis_audit.csv")
)


# ============================================================
# 15. REPORT TESTUALE
# ============================================================

report_path <- file.path(
  OUTPUT_REPORTS, "integrated_gene_module_analysis_report.txt"
)

sink(report_path)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("INTEGRATED GENE MODULE ANALYSIS - SCRIPT 33\n")
cat("============================================================\n\n")

cat("PROJECT:", PROJECT_ROOT, "\n")
cat("FDR CUTOFF:", FDR_CUTOFF, "\n\n")

cat("INPUT FILES\n")
cat("------------------------------------------------------------\n")
cat("DE 1:", de_file_1, "\n")
cat("DE 2:", de_file_2, "\n")
cat("SIGNATURES:", signature_path, "\n\n")

cat("GLOBAL DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")
print(input_audit, row.names = FALSE)

cat("\nMODULE SUMMARY\n")
cat("------------------------------------------------------------\n")
print(module_summary, row.names = FALSE)

cat("\nMODULE GENES WITH SIGNIFICANT RESULTS\n")
cat("------------------------------------------------------------\n")

for (module_name in names(GENE_MODULES)) {
  cat("\nMODULE:", module_name, "\n")

  selected <- module_comparison[
    module_comparison$module == module_name &
      module_comparison$significant_either,
    ,
    drop = FALSE
  ]

  if (nrow(selected) == 0) {
    cat("Nessun gene significativo rilevato.\n")
  } else {
    selected$minimum_FDR_report <- pmin(
      selected$FDR_siNRP1_3,
      selected$FDR_siNRP1_5,
      na.rm = TRUE
    )

    selected <- selected[
      order(selected$minimum_FDR_report, na.last = TRUE),
      ,
      drop = FALSE
    ]

    print(
      selected[, c(
        "gene_symbol",
        "logFC_siNRP1_3",
        "FDR_siNRP1_3",
        "logFC_siNRP1_5",
        "FDR_siNRP1_5",
        "direction_concordance",
        "joint_class"
      )],
      row.names = FALSE
    )
  }
}

cat("\nSIGNATURE-MODULE OVERLAP\n")
cat("------------------------------------------------------------\n")
if (nrow(signature_module_summary) == 0) {
  cat("Nessuna sovrapposizione interpretabile disponibile.\n")
} else {
  print(signature_module_summary, row.names = FALSE)
}

cat("\nFINAL AUDIT\n")
cat("------------------------------------------------------------\n")
print(final_audit, row.names = FALSE)

cat("\nINTERPRETATION NOTES\n")
cat("------------------------------------------------------------\n")
cat("1. I moduli sono liste candidate definite a priori.\n")
cat("2. L'appartenenza a un modulo non dimostra attivazione di pathway.\n")
cat("3. Il logFC descrive l'espressione relativa rispetto a siNT.\n")
cat("4. La concordanza indica lo stesso segno del logFC nei due contrasti.\n")
cat("5. I simboli duplicati sono risolti scegliendo il record con FDR minore.\n")
cat("6. I record duplicati originali sono documentati nell'audit.\n")
cat("7. Le analisi sono descrittive e non dimostrano causalita'.\n")
cat("8. Le liste candidate non costituiscono un test formale di gene set.\n")
cat("9. Le firme enrichment sono integrate solo se contengono simboli leggibili.\n")

cat("\n============================================================\n")
cat("END OF REPORT\n")
cat("============================================================\n")

sink()


# ============================================================
# 16. SESSION INFO
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_AUDIT, "sessionInfo.txt")
)


# ============================================================
# 17. COMPLETAMENTO
# ============================================================

log_message("SCRIPT 33 COMPLETATO.")
log_message("Output directory:", OUTPUT_DIR)
log_message("Moduli analizzati:", length(GENE_MODULES))
log_message("Geni candidati unici:", length(all_module_genes))
log_message("Record integrati:", nrow(module_comparison))
log_message(
  "Geni significativi in entrambi:",
  sum(module_comparison$significant_both, na.rm = TRUE)
)
log_message(
  "Geni concordanti tra i significativi in entrambi:",
  sum(
    module_comparison$significant_both &
      module_comparison$direction_concordance == "CONCORDANT",
    na.rm = TRUE
  )
)
log_message(
  "Geni discordanti tra i significativi in entrambi:",
  sum(
    module_comparison$significant_both &
      module_comparison$direction_concordance == "DISCORDANT",
    na.rm = TRUE
  )
)

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - SCRIPT 33 FINISHED\n")
cat("============================================================\n")
cat("Output:", OUTPUT_DIR, "\n")
