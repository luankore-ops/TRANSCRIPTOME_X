
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 41 - INTEGRATED CANDIDATE GENE VALIDATION
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 41 - INTEGRATED CANDIDATE GENE VALIDATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(project_dir, "results")

input_files <- c(
  network = file.path(
    results_dir,
    "network_candidate_prioritization/tables/candidate_gene_priority_ranking.csv"
  ),
  de3 = file.path(
    results_dir,
    "normalization_sensitivity/DE_TMM_siNRP1_3_vs_siNT.csv"
  ),
  de5 = file.path(
    results_dir,
    "normalization_sensitivity/DE_TMM_siNRP1_5_vs_siNT.csv"
  ),
  pathway = file.path(
    results_dir,
    "pathway_concordance_specificity/tables/pathway_gene_evidence.csv"
  ),
  module = file.path(
    results_dir,
    "gene_pathway_module_network/tables/gene_module_edges.csv"
  )
)

output_dir <- file.path(
  results_dir,
  "integrated_candidate_gene_validation"
)

table_dir <- file.path(output_dir, "tables")
report_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (d in c(table_dir, report_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. FUNCTIONS
# ------------------------------------------------------------

read_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN", "N/A", "#N/A")
  )
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

safe_numeric <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "N/A", "#N/A")] <- NA_character_
  
  # Decimal comma support only when no decimal point is present.
  comma_decimal <- !is.na(x) &
    grepl(",", x, fixed = TRUE) &
    !grepl(".", x, fixed = TRUE)
  
  x[comma_decimal] <- gsub(
    ",", ".", x[comma_decimal], fixed = TRUE
  )
  
  suppressWarnings(as.numeric(x))
}

find_column <- function(data, candidates, label) {
  
  normalized <- tolower(
    gsub("[^a-z0-9]", "", names(data))
  )
  
  targets <- tolower(
    gsub("[^a-z0-9]", "", candidates)
  )
  
  idx <- match(targets, normalized)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    stop(
      "Missing column in ", label,
      ". Available: ", paste(names(data), collapse = ", ")
    )
  }
  
  names(data)[idx[1]]
}

write_csv <- function(data, path) {
  write.csv(
    data,
    path,
    row.names = FALSE,
    na = ""
  )
}

rank_score <- function(x) {
  
  x <- safe_numeric(x)
  out <- rep(0, length(x))
  valid <- is.finite(x)
  
  if (sum(valid) == 1) {
    out[valid] <- 1
  } else if (sum(valid) > 1) {
    r <- rank(x[valid], ties.method = "average")
    out[valid] <- (r - 1) / (length(r) - 1)
  }
  
  out
}

# ------------------------------------------------------------
# 3. LOAD INPUTS
# ------------------------------------------------------------

cat("Loading existing analysis results...\n")

network_raw <- read_csv(input_files["network"])
de3_raw <- read_csv(input_files["de3"])
de5_raw <- read_csv(input_files["de5"])
pathway_raw <- read_csv(input_files["pathway"])
module_raw <- read_csv(input_files["module"])

cat("All input files loaded.\n\n")

cat("INPUT DIMENSIONS\n")
cat("Network rows:", nrow(network_raw), "\n")
cat("DE siNRP1.3 rows:", nrow(de3_raw), "\n")
cat("DE siNRP1.5 rows:", nrow(de5_raw), "\n")
cat("Pathway rows:", nrow(pathway_raw), "\n")
cat("Module rows:", nrow(module_raw), "\n\n")

# ------------------------------------------------------------
# 4. NETWORK RANKING
# ------------------------------------------------------------

network_gene_col <- find_column(
  network_raw,
  c("gene_symbol", "gene", "symbol"),
  "network ranking"
)

network_score_col <- find_column(
  network_raw,
  c("priority_score"),
  "network ranking"
)

network_class_col <- find_column(
  network_raw,
  c("priority_class"),
  "network ranking"
)

network <- data.frame(
  gene_symbol = clean_gene(network_raw[[network_gene_col]]),
  network_score = safe_numeric(network_raw[[network_score_col]]),
  network_class = as.character(network_raw[[network_class_col]]),
  stringsAsFactors = FALSE
)

network <- network[
  !is.na(network$gene_symbol),
  ,
  drop = FALSE
]

network <- network[
  !duplicated(network$gene_symbol),
  ,
  drop = FALSE
]

if (any(
  !is.na(network$network_score) &
  (
    !is.finite(network$network_score) |
    network$network_score < 0 |
    network$network_score > 100
  )
)) {
  stop("Invalid network priority scores.")
}

cat("Unique network candidates:", nrow(network), "\n\n")

# ------------------------------------------------------------
# 5. DIFFERENTIAL EXPRESSION
# ------------------------------------------------------------

prepare_de <- function(data, comparison_name) {
  
  gene_col <- find_column(
    data,
    c(
      "Ensembl.114.Gene.Name",
      "gene_symbol",
      "symbol",
      "gene_name"
    ),
    comparison_name
  )
  
  logfc_col <- find_column(
    data,
    c("logFC", "log2FoldChange"),
    comparison_name
  )
  
  # IMPORTANT:
  # Use the exact FDR column, not a flexible column search.
  if (!"FDR" %in% names(data)) {
    stop(
      "Exact FDR column not found in ", comparison_name,
      ". Available columns: ",
      paste(names(data), collapse = ", ")
    )
  }
  
  raw_fdr <- data[["FDR"]]
  
  result <- data.frame(
    gene_symbol = clean_gene(data[[gene_col]]),
    logFC = safe_numeric(data[[logfc_col]]),
    FDR = safe_numeric(raw_fdr),
    stringsAsFactors = FALSE
  )
  
  # Keep only rows with a usable gene symbol.
  result <- result[
    !is.na(result$gene_symbol),
    ,
    drop = FALSE
  ]
  
  # Missing FDR is allowed.
  # Non-missing values outside [0,1] are invalid.
  invalid_fdr <- !is.na(result$FDR) & (
    !is.finite(result$FDR) |
      result$FDR < 0 |
      result$FDR > 1
  )
  
  # Values that could not be converted to numeric are
  # treated as missing, but explicitly counted and audited.
  raw_fdr_after_gene_filter <- raw_fdr[
    !is.na(clean_gene(data[[gene_col]]))
  ]
  
  converted_fdr <- safe_numeric(raw_fdr_after_gene_filter)
  
  conversion_failures <- sum(
    !is.na(raw_fdr_after_gene_filter) &
      trimws(as.character(raw_fdr_after_gene_filter)) != "" &
      is.na(converted_fdr)
  )
  
  cat("\nFDR DIAGNOSTICS -", comparison_name, "\n")
  cat("FDR column: FDR\n")
  cat("Rows with gene symbol:", nrow(result), "\n")
  cat("Missing FDR:", sum(is.na(result$FDR)), "\n")
  cat("Non-finite FDR:", sum(
    !is.na(result$FDR) & !is.finite(result$FDR)
  ), "\n")
  cat("FDR below zero:", sum(
    !is.na(result$FDR) & result$FDR < 0
  ), "\n")
  cat("FDR above one:", sum(
    !is.na(result$FDR) & result$FDR > 1
  ), "\n")
  cat("Unparseable FDR entries:", conversion_failures, "\n")
  
  if (any(invalid_fdr)) {
    bad <- result[invalid_fdr, c("gene_symbol", "FDR")]
    write_csv(
      bad,
      file.path(
        audit_dir,
        paste0("invalid_FDR_", comparison_name, ".csv")
      )
    )
    
    stop(
      "Invalid numeric FDR values in ", comparison_name,
      ". Invalid rows: ", sum(invalid_fdr),
      ". See audit/invalid_FDR_", comparison_name, ".csv"
    )
  }
  
  invalid_logfc <- !is.na(result$logFC) &
    !is.finite(result$logFC)
  
  if (any(invalid_logfc)) {
    stop(
      "Invalid logFC values in ", comparison_name,
      ": ", sum(invalid_logfc)
    )
  }
  
  result$.missing_fdr <- is.na(result$FDR)
  
  result <- result[
    order(
      result$.missing_fdr,
      result$FDR,
      -abs(result$logFC),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  result <- result[
    !duplicated(result$gene_symbol),
    ,
    drop = FALSE
  ]
  
  result$.missing_fdr <- NULL
  
  names(result)[names(result) == "logFC"] <-
    paste0("logFC_", comparison_name)
  
  names(result)[names(result) == "FDR"] <-
    paste0("FDR_", comparison_name)
  
  result
}

cat("Processing differential expression...\n")

de3 <- prepare_de(de3_raw, "siNRP1_3")
de5 <- prepare_de(de5_raw, "siNRP1_5")

cat("\nDE processing completed.\n")
cat("Unique genes siNRP1.3:", nrow(de3), "\n")
cat("Unique genes siNRP1.5:", nrow(de5), "\n\n")

# ------------------------------------------------------------
# 6. NRP1 INTEGRITY CHECK
# ------------------------------------------------------------

check_nrp1 <- function(data, comparison, expected_logFC, expected_FDR) {
  
  idx <- which(data$gene_symbol == "NRP1")
  
  if (length(idx) != 1) {
    stop("NRP1 missing or duplicated in ", comparison)
  }
  
  logfc_col <- paste0("logFC_", comparison)
  fdr_col <- paste0("FDR_", comparison)
  
  if (abs(data[[logfc_col]][idx] - expected_logFC) > 1e-6) {
    stop("NRP1 logFC integrity check failed: ", comparison)
  }
  
  if (abs(data[[fdr_col]][idx] - expected_FDR) > 1e-6) {
    stop("NRP1 FDR integrity check failed: ", comparison)
  }
  
  invisible(TRUE)
}

check_nrp1(
  de3,
  "siNRP1_3",
  -0.692876967820272,
  0.0247835265929449
)

check_nrp1(
  de5,
  "siNRP1_5",
  -4.87695296739976,
  3.02363667601604e-07
)

cat("NRP1 integrity checks passed.\n\n")

# ------------------------------------------------------------
# 7. PATHWAY EVIDENCE
# ------------------------------------------------------------

if (!all(c("gene_symbol", "pathway_id") %in% names(pathway_raw))) {
  stop(
    "Pathway input must contain gene_symbol and pathway_id. Found: ",
    paste(names(pathway_raw), collapse = ", ")
  )
}

pathway <- unique(data.frame(
  gene_symbol = clean_gene(pathway_raw$gene_symbol),
  pathway_id = trimws(as.character(pathway_raw$pathway_id)),
  stringsAsFactors = FALSE
))

pathway <- pathway[
  !is.na(pathway$gene_symbol) &
    !is.na(pathway$pathway_id) &
    pathway$pathway_id != "",
  ,
  drop = FALSE
]

pathway_split <- split(
  pathway$pathway_id,
  pathway$gene_symbol
)

pathway_summary <- data.frame(
  gene_symbol = names(pathway_split),
  pathway_count = vapply(
    pathway_split,
    function(x) length(unique(x)),
    integer(1)
  ),
  pathway_list = vapply(
    pathway_split,
    function(x) paste(sort(unique(x)), collapse = "; "),
    character(1)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 8. MODULE EVIDENCE
# ------------------------------------------------------------

module_gene_col <- find_column(
  module_raw,
  c("gene_symbol", "gene", "symbol"),
  "gene module edges"
)

module_name_col <- find_column(
  module_raw,
  c("module", "module_name"),
  "gene module edges"
)

module <- unique(data.frame(
  gene_symbol = clean_gene(module_raw[[module_gene_col]]),
  module = trimws(as.character(module_raw[[module_name_col]])),
  stringsAsFactors = FALSE
))

module <- module[
  !is.na(module$gene_symbol) &
    !is.na(module$module) &
    module$module != "" &
    module$module != "NO_DEFINED_MODULE",
  ,
  drop = FALSE
]

module_split <- split(
  module$module,
  module$gene_symbol
)

module_summary <- data.frame(
  gene_symbol = names(module_split),
  module_count = vapply(
    module_split,
    function(x) length(unique(x)),
    integer(1)
  ),
  module_list = vapply(
    module_split,
    function(x) paste(sort(unique(x)), collapse = "; "),
    character(1)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9. INTEGRATE EVIDENCE
# ------------------------------------------------------------

cat("Integrating evidence layers...\n")

integrated <- merge(
  network, de3,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

integrated <- merge(
  integrated, de5,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

integrated <- merge(
  integrated, pathway_summary,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

integrated <- merge(
  integrated, module_summary,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

integrated$pathway_count[
  is.na(integrated$pathway_count)
] <- 0

integrated$module_count[
  is.na(integrated$module_count)
] <- 0

integrated$pathway_list[
  is.na(integrated$pathway_list)
] <- ""

integrated$module_list[
  is.na(integrated$module_list)
] <- ""

if (anyDuplicated(integrated$gene_symbol)) {
  stop("Duplicated genes after evidence integration.")
}

if (nrow(integrated) != nrow(network)) {
  stop("Integrated candidate count differs from network ranking.")
}

cat("Integrated candidates:", nrow(integrated), "\n\n")

# ------------------------------------------------------------
# 10. DIFFERENTIAL EXPRESSION CLASSIFICATION
# ------------------------------------------------------------

integrated$significant_siNRP1_3 <- (
  !is.na(integrated$FDR_siNRP1_3) &
    integrated$FDR_siNRP1_3 < 0.05
)

integrated$significant_siNRP1_5 <- (
  !is.na(integrated$FDR_siNRP1_5) &
    integrated$FDR_siNRP1_5 < 0.05
)

integrated$significant_both <- (
  integrated$significant_siNRP1_3 &
    integrated$significant_siNRP1_5
)

integrated$significant_either <- (
  integrated$significant_siNRP1_3 |
    integrated$significant_siNRP1_5
)

integrated$direction_siNRP1_3 <- ifelse(
  is.na(integrated$logFC_siNRP1_3),
  "NOT_AVAILABLE",
  ifelse(
    integrated$logFC_siNRP1_3 > 0, "UP",
    ifelse(integrated$logFC_siNRP1_3 < 0, "DOWN", "ZERO")
  )
)

integrated$direction_siNRP1_5 <- ifelse(
  is.na(integrated$logFC_siNRP1_5),
  "NOT_AVAILABLE",
  ifelse(
    integrated$logFC_siNRP1_5 > 0, "UP",
    ifelse(integrated$logFC_siNRP1_5 < 0, "DOWN", "ZERO")
  )
)

integrated$direction_concordance <- ifelse(
  integrated$significant_both &
    integrated$direction_siNRP1_3 ==
    integrated$direction_siNRP1_5,
  "CONCORDANT",
  ifelse(
    integrated$significant_both,
    "DISCORDANT",
    ifelse(
      integrated$significant_either,
      "SINGLE_CONTRAST_SUPPORT",
      "NOT_SIGNIFICANT"
    )
  )
)

# ------------------------------------------------------------
# 11. EVIDENCE LAYERS
# ------------------------------------------------------------

integrated$has_network_evidence <- !is.na(
  integrated$network_score
)

integrated$has_DE_evidence <- integrated$significant_either

integrated$has_pathway_evidence <- (
  integrated$pathway_count > 0
)

integrated$has_module_evidence <- (
  integrated$module_count > 0
)

integrated$evidence_layer_count <- (
  as.integer(integrated$has_network_evidence) +
    as.integer(integrated$has_DE_evidence) +
    as.integer(integrated$has_pathway_evidence) +
    as.integer(integrated$has_module_evidence)
)

# ------------------------------------------------------------
# 12. COMPONENT SCORES
# ------------------------------------------------------------

integrated$network_component <- (
  pmax(0, pmin(100, integrated$network_score)) / 100
)

integrated$DE_component <- ifelse(
  integrated$significant_both,
  1,
  ifelse(integrated$significant_either, 0.5, 0)
)

integrated$concordance_component <- ifelse(
  integrated$direction_concordance == "CONCORDANT",
  1,
  ifelse(
    integrated$direction_concordance == "SINGLE_CONTRAST_SUPPORT",
    0.5,
    0
  )
)

integrated$pathway_component <- rank_score(
  integrated$pathway_count
)

integrated$module_component <- rank_score(
  integrated$module_count
)

# ------------------------------------------------------------
# 13. INTEGRATED SCORE
# ------------------------------------------------------------

integrated$validation_score <- 100 * (
  0.25 * integrated$network_component +
    0.30 * integrated$DE_component +
    0.20 * integrated$concordance_component +
    0.15 * integrated$pathway_component +
    0.10 * integrated$module_component
)

# ------------------------------------------------------------
# 14. VALIDATION CLASS
# ------------------------------------------------------------

integrated$validation_class <- ifelse(
  integrated$significant_both &
    integrated$direction_concordance == "CONCORDANT" &
    integrated$pathway_count > 0 &
    integrated$module_count > 0,
  "ROBUST_CROSS_LAYER",
  ifelse(
    integrated$significant_either &
      integrated$pathway_count > 0 &
      integrated$module_count > 0,
    "MULTI_LAYER_SUPPORTED",
    ifelse(
      integrated$significant_either,
      "TRANSCRIPTOMICALLY_SUPPORTED",
      ifelse(
        integrated$network_class == "HIGH",
        "NETWORK_PRIORITY_ONLY",
        "EXPLORATORY"
      )
    )
  )
)

# ------------------------------------------------------------
# 15. RANKING
# ------------------------------------------------------------

integrated <- integrated[
  order(
    -integrated$validation_score,
    -as.integer(integrated$significant_both),
    -integrated$pathway_count,
    -integrated$module_count,
    integrated$gene_symbol
  ),
  ,
  drop = FALSE
]

rownames(integrated) <- NULL
integrated$validation_rank <- seq_len(nrow(integrated))

# ------------------------------------------------------------
# 16. SHORTLISTS
# ------------------------------------------------------------

shortlist <- integrated[
  integrated$validation_class %in% c(
    "ROBUST_CROSS_LAYER",
    "MULTI_LAYER_SUPPORTED"
  ),
  ,
  drop = FALSE
]

robust_candidates <- integrated[
  integrated$validation_class == "ROBUST_CROSS_LAYER",
  ,
  drop = FALSE
]

top_50 <- head(integrated, 50)

focused_genes <- c(
  "THBS1", "CCN3", "TGFB2", "BMP4", "CAV1",
  "NRP1", "LYN", "TGFBR3", "ENG", "TGFBR2", "CCN2"
)

focused_panel <- integrated[
  integrated$gene_symbol %in% focused_genes,
  ,
  drop = FALSE
]

focused_panel <- focused_panel[
  order(-focused_panel$validation_score),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 17. EXPORT TABLES
# ------------------------------------------------------------

cat("Exporting results...\n")

write_csv(
  integrated,
  file.path(table_dir, "integrated_candidate_gene_validation.csv")
)

write_csv(
  shortlist,
  file.path(table_dir, "final_candidate_gene_shortlist.csv")
)

write_csv(
  robust_candidates,
  file.path(table_dir, "robust_cross_layer_candidates.csv")
)

write_csv(
  top_50,
  file.path(table_dir, "top_50_integrated_candidates.csv")
)

write_csv(
  focused_panel,
  file.path(table_dir, "focused_candidate_gene_panel.csv")
)

# ------------------------------------------------------------
# 18. CLASS SUMMARY
# ------------------------------------------------------------

class_order <- c(
  "ROBUST_CROSS_LAYER",
  "MULTI_LAYER_SUPPORTED",
  "TRANSCRIPTOMICALLY_SUPPORTED",
  "NETWORK_PRIORITY_ONLY",
  "EXPLORATORY"
)

class_summary <- as.data.frame(
  table(
    factor(
      integrated$validation_class,
      levels = class_order
    )
  ),
  stringsAsFactors = FALSE
)

names(class_summary) <- c(
  "validation_class",
  "gene_count"
)

class_summary$percentage <- (
  class_summary$gene_count / nrow(integrated)
) * 100

write_csv(
  class_summary,
  file.path(table_dir, "integrated_validation_class_summary.csv")
)

# ------------------------------------------------------------
# 19. EVIDENCE SUMMARY
# ------------------------------------------------------------

evidence_summary <- data.frame(
  evidence_type = c(
    "Network ranking available",
    "Significant in siNRP1.3",
    "Significant in siNRP1.5",
    "Significant in both",
    "Concordant in both",
    "Discordant in both",
    "Pathway-supported",
    "Module-supported",
    "All four evidence layers"
  ),
  gene_count = c(
    sum(integrated$has_network_evidence),
    sum(integrated$significant_siNRP1_3),
    sum(integrated$significant_siNRP1_5),
    sum(integrated$significant_both),
    sum(integrated$direction_concordance == "CONCORDANT"),
    sum(integrated$direction_concordance == "DISCORDANT"),
    sum(integrated$has_pathway_evidence),
    sum(integrated$has_module_evidence),
    sum(integrated$evidence_layer_count == 4)
  )
)

write_csv(
  evidence_summary,
  file.path(table_dir, "integrated_evidence_summary.csv")
)

# ------------------------------------------------------------
# 20. TOP CANDIDATES REPORT
# ------------------------------------------------------------

top10 <- head(integrated, 10)

top10_lines <- sprintf(
  "%d. %s | Score %.2f | %s | %s | Pathways %d | Modules %d",
  top10$validation_rank,
  top10$gene_symbol,
  top10$validation_score,
  top10$validation_class,
  top10$direction_concordance,
  top10$pathway_count,
  top10$module_count
)

report <- c(
  "TRANSCRIPTOME-X",
  "SCRIPT 41 - INTEGRATED CANDIDATE GENE VALIDATION",
  "============================================================",
  "",
  paste("Integrated candidates:", nrow(integrated)),
  paste("Final shortlist:", nrow(shortlist)),
  paste("Robust cross-layer candidates:", nrow(robust_candidates)),
  "",
  "CLASS DISTRIBUTION",
  "------------------",
  paste(
    class_summary$validation_class,
    class_summary$gene_count,
    sep = ": "
  ),
  "",
  "TOP 10 CANDIDATES",
  "-----------------",
  top10_lines,
  "",
  "SCORING FRAMEWORK",
  "-----------------",
  "Network evidence: 25%",
  "Differential expression: 30%",
  "Directional agreement: 20%",
  "Pathway support: 15%",
  "Module support: 10%",
  "",
  "INTERPRETATION",
  "--------------",
  "The integrated score is a descriptive prioritization score.",
  "It is not a statistical significance measure or causal estimate.",
  "Pathway and module annotations are overlapping evidence layers.",
  "Missing FDR values are not considered statistically significant.",
  "Concordance across knockdown conditions does not establish causality."
)

writeLines(
  report,
  file.path(report_dir, "integrated_candidate_gene_validation_report.txt")
)

# ------------------------------------------------------------
# 21. AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "network_input_rows",
    "DE_siNRP1_3_input_rows",
    "DE_siNRP1_5_input_rows",
    "DE_siNRP1_3_unique_genes",
    "DE_siNRP1_5_unique_genes",
    "DE_siNRP1_3_missing_FDR",
    "DE_siNRP1_5_missing_FDR",
    "DE_siNRP1_3_missing_logFC",
    "DE_siNRP1_5_missing_logFC",
    "pathway_input_rows",
    "unique_pathway_gene_pairs",
    "module_input_rows",
    "unique_module_gene_pairs",
    "integrated_candidates",
    "shortlist_candidates",
    "robust_candidates",
    "significant_both",
    "concordant_both",
    "discordant_both"
  ),
  value = c(
    nrow(network_raw),
    nrow(de3_raw),
    nrow(de5_raw),
    nrow(de3),
    nrow(de5),
    sum(is.na(de3$FDR_siNRP1_3)),
    sum(is.na(de5$FDR_siNRP1_5)),
    sum(is.na(de3$logFC_siNRP1_3)),
    sum(is.na(de5$logFC_siNRP1_5)),
    nrow(pathway_raw),
    nrow(pathway),
    nrow(module_raw),
    nrow(module),
    nrow(integrated),
    nrow(shortlist),
    nrow(robust_candidates),
    sum(integrated$significant_both),
    sum(integrated$direction_concordance == "CONCORDANT"),
    sum(integrated$direction_concordance == "DISCORDANT")
  )
)

write_csv(
  audit,
  file.path(audit_dir, "script41_audit.csv")
)

# ------------------------------------------------------------
# 22. FINAL VALIDATION
# ------------------------------------------------------------

if (anyDuplicated(integrated$gene_symbol)) {
  stop("Duplicate gene symbols in final output.")
}

if (anyDuplicated(integrated$validation_rank)) {
  stop("Duplicate validation ranks.")
}

if (any(!is.finite(integrated$validation_score))) {
  stop("Non-finite validation scores.")
}

if (any(
  integrated$validation_score < 0 |
  integrated$validation_score > 100
)) {
  stop("Validation scores outside 0-100.")
}

for (col in c("FDR_siNRP1_3", "FDR_siNRP1_5")) {
  
  x <- integrated[[col]]
  
  invalid <- !is.na(x) & (
    !is.finite(x) | x < 0 | x > 1
  )
  
  if (any(invalid)) {
    stop("Invalid final FDR values in ", col)
  }
}

if (any(
  is.na(integrated$FDR_siNRP1_3) &
  integrated$significant_siNRP1_3
)) {
  stop("Missing FDR incorrectly classified as significant.")
}

if (any(
  is.na(integrated$FDR_siNRP1_5) &
  integrated$significant_siNRP1_5
)) {
  stop("Missing FDR incorrectly classified as significant.")
}

# ------------------------------------------------------------
# 23. COMPLETION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 41 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Integrated candidates:", nrow(integrated), "\n")
cat("Final shortlist:", nrow(shortlist), "\n")
cat("Robust cross-layer candidates:", nrow(robust_candidates), "\n\n")

cat("CLASS DISTRIBUTION\n")
print(class_summary, row.names = FALSE)

cat("\nTOP 10 CANDIDATES\n")
print(
  top10[
    ,
    c(
      "validation_rank",
      "gene_symbol",
      "validation_score",
      "validation_class",
      "direction_concordance"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nFOCUSED CANDIDATE PANEL\n")
print(
  focused_panel[
    ,
    c(
      "gene_symbol",
      "validation_score",
      "validation_class",
      "logFC_siNRP1_3",
      "FDR_siNRP1_3",
      "logFC_siNRP1_5",
      "FDR_siNRP1_5"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
