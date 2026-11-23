
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 42 - FINAL BIOLOGICAL INTERPRETATION
# AND CANDIDATE PRIORITIZATION
#
# INPUT:
# Script 41 - Integrated candidate gene validation
#
# OBJECTIVES:
# 1. Consolidate the 1,913 integrated candidates.
# 2. Prioritize the 79 multi-layer supported candidates.
# 3. Characterize the 23 robust cross-layer candidates.
# 4. Compare expression patterns across both NRP1 knockdowns.
# 5. Summarize pathway and module evidence.
# 6. Generate final candidate tables and interpretation reports.
#
# No new differential expression or enrichment analysis.
# No causal or mechanistic claims.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 42 - FINAL BIOLOGICAL INTERPRETATION\n")
cat("AND CANDIDATE PRIORITIZATION\n")
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

input_file <- file.path(
  results_dir,
  "integrated_candidate_gene_validation",
  "tables",
  "integrated_candidate_gene_validation.csv"
)

output_dir <- file.path(
  results_dir,
  "final_biological_interpretation"
)

table_dir <- file.path(output_dir, "tables")
report_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (d in c(table_dir, report_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) {
  stop("Script 41 integrated candidate table not found:\n", input_file)
}


# ------------------------------------------------------------
# 2. FUNCTIONS
# ------------------------------------------------------------

read_csv <- function(path) {
  
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN", "N/A")
  )
}

write_csv <- function(data, path) {
  
  write.csv(
    data,
    path,
    row.names = FALSE,
    na = ""
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

find_column <- function(data, candidates, required = TRUE) {
  
  normalized <- tolower(
    gsub("[^a-z0-9]", "", names(data))
  )
  
  targets <- tolower(
    gsub("[^a-z0-9]", "", candidates)
  )
  
  idx <- match(targets, normalized)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    
    if (required) {
      stop(
        "Required column not found. Candidates: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns: ",
        paste(names(data), collapse = ", ")
      )
    }
    
    return(NULL)
  }
  
  names(data)[idx[1]]
}

assert_columns <- function(data, columns, label) {
  
  missing <- setdiff(columns, names(data))
  
  if (length(missing) > 0) {
    stop(
      "Missing columns in ", label, ": ",
      paste(missing, collapse = ", ")
    )
  }
}

split_annotations <- function(x) {
  
  if (is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(character(0))
  }
  
  values <- trimws(
    unlist(
      strsplit(as.character(x), ";", fixed = TRUE)
    )
  )
  
  values <- values[
    !is.na(values) & values != ""
  ]
  
  unique(values)
}

annotation_frequency <- function(data, annotation_column) {
  
  if (nrow(data) == 0) {
    return(data.frame(
      annotation = character(0),
      candidate_count = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  annotations <- lapply(
    data[[annotation_column]],
    split_annotations
  )
  
  all_annotations <- unlist(
    annotations,
    use.names = FALSE
  )
  
  if (length(all_annotations) == 0) {
    return(data.frame(
      annotation = character(0),
      candidate_count = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  counts <- sort(
    table(all_annotations),
    decreasing = TRUE
  )
  
  result <- data.frame(
    annotation = names(counts),
    candidate_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
  
  result$percentage <- (
    result$candidate_count / nrow(data)
  ) * 100
  
  result
}

format_gene_list <- function(x) {
  
  if (length(x) == 0) {
    return("None")
  }
  
  paste(x, collapse = ", ")
}


# ------------------------------------------------------------
# 3. LOAD SCRIPT 41 RESULTS
# ------------------------------------------------------------

cat("Loading Script 41 integrated results...\n")

integrated <- read_csv(input_file)

cat("Input loaded.\n")
cat("Rows:", nrow(integrated), "\n")
cat("Columns:", ncol(integrated), "\n\n")


# ------------------------------------------------------------
# 4. INPUT VALIDATION
# ------------------------------------------------------------

required_columns <- c(
  "gene_symbol",
  "network_score",
  "network_class",
  "logFC_siNRP1_3",
  "FDR_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_5",
  "significant_siNRP1_3",
  "significant_siNRP1_5",
  "significant_both",
  "direction_concordance",
  "pathway_count",
  "pathway_list",
  "module_count",
  "module_list",
  "validation_score",
  "validation_class",
  "validation_rank"
)

assert_columns(
  integrated,
  required_columns,
  "Script 41 integrated results"
)

integrated$gene_symbol <- trimws(
  as.character(integrated$gene_symbol)
)

if (any(
  is.na(integrated$gene_symbol) |
  integrated$gene_symbol == ""
)) {
  stop("Missing gene symbols detected.")
}

if (anyDuplicated(integrated$gene_symbol)) {
  stop("Duplicated gene symbols detected.")
}

numeric_columns <- c(
  "network_score",
  "logFC_siNRP1_3",
  "FDR_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_5",
  "pathway_count",
  "module_count",
  "validation_score",
  "validation_rank"
)

for (col in numeric_columns) {
  integrated[[col]] <- safe_numeric(integrated[[col]])
}

for (col in c("FDR_siNRP1_3", "FDR_siNRP1_5")) {
  
  x <- integrated[[col]]
  
  invalid <- !is.na(x) & (
    !is.finite(x) | x < 0 | x > 1
  )
  
  if (any(invalid)) {
    stop("Invalid FDR values detected in ", col)
  }
}

if (any(
  !is.finite(integrated$validation_score) |
  integrated$validation_score < 0 |
  integrated$validation_score > 100
)) {
  stop("Invalid integrated validation scores.")
}

cat("Input integrity checks passed.\n\n")


# ------------------------------------------------------------
# 5. RECONSTRUCT FINAL CANDIDATE GROUPS
# ------------------------------------------------------------

cat("Reconstructing candidate groups...\n")

robust <- integrated[
  integrated$validation_class == "ROBUST_CROSS_LAYER",
  ,
  drop = FALSE
]

multi_layer <- integrated[
  integrated$validation_class == "MULTI_LAYER_SUPPORTED",
  ,
  drop = FALSE
]

transcriptomic <- integrated[
  integrated$validation_class ==
    "TRANSCRIPTOMICALLY_SUPPORTED",
  ,
  drop = FALSE
]

network_only <- integrated[
  integrated$validation_class ==
    "NETWORK_PRIORITY_ONLY",
  ,
  drop = FALSE
]

exploratory <- integrated[
  integrated$validation_class == "EXPLORATORY",
  ,
  drop = FALSE
]

shortlist <- integrated[
  integrated$validation_class %in% c(
    "ROBUST_CROSS_LAYER",
    "MULTI_LAYER_SUPPORTED"
  ),
  ,
  drop = FALSE
]

if (nrow(integrated) != 1913) {
  warning(
    "The number of integrated candidates differs from the ",
    "previous Script 41 execution."
  )
}

if (nrow(shortlist) != 79) {
  warning(
    "The shortlist size differs from the previous Script 41 execution."
  )
}

cat("Integrated candidates:", nrow(integrated), "\n")
cat("Robust cross-layer:", nrow(robust), "\n")
cat("Multi-layer supported:", nrow(multi_layer), "\n")
cat("Final shortlist:", nrow(shortlist), "\n\n")


# ------------------------------------------------------------
# 6. EXPRESSION PATTERN CLASSIFICATION
# ------------------------------------------------------------

cat("Characterizing expression patterns...\n")

integrated$expression_pattern <- "NOT_SIGNIFICANT"

both_significant <- (
  integrated$significant_siNRP1_3 &
    integrated$significant_siNRP1_5
)

only_3 <- (
  integrated$significant_siNRP1_3 &
    !integrated$significant_siNRP1_5
)

only_5 <- (
  !integrated$significant_siNRP1_3 &
    integrated$significant_siNRP1_5
)

integrated$expression_pattern[only_3] <-
  "SIGNIFICANT_ONLY_siNRP1_3"

integrated$expression_pattern[only_5] <-
  "SIGNIFICANT_ONLY_siNRP1_5"

integrated$expression_pattern[
  both_significant &
    integrated$logFC_siNRP1_3 < 0 &
    integrated$logFC_siNRP1_5 < 0
] <- "DOWN_BOTH"

integrated$expression_pattern[
  both_significant &
    integrated$logFC_siNRP1_3 > 0 &
    integrated$logFC_siNRP1_5 > 0
] <- "UP_BOTH"

integrated$expression_pattern[
  both_significant &
    integrated$logFC_siNRP1_3 *
    integrated$logFC_siNRP1_5 < 0
] <- "OPPOSITE_DIRECTIONS"

integrated$expression_pattern[
  both_significant &
    (
      integrated$logFC_siNRP1_3 == 0 |
        integrated$logFC_siNRP1_5 == 0
    )
] <- "ZERO_LOGFC_IN_ONE_CONTRAST"


# ------------------------------------------------------------
# 7. FINAL CANDIDATE PRIORITIZATION
# ------------------------------------------------------------

integrated$final_priority <- ifelse(
  integrated$validation_class == "ROBUST_CROSS_LAYER",
  "PRIORITY_1_ROBUST",
  ifelse(
    integrated$validation_class == "MULTI_LAYER_SUPPORTED",
    "PRIORITY_2_MULTI_LAYER",
    ifelse(
      integrated$validation_class ==
        "TRANSCRIPTOMICALLY_SUPPORTED",
      "PRIORITY_3_TRANSCRIPTOMIC",
      ifelse(
        integrated$validation_class ==
          "NETWORK_PRIORITY_ONLY",
        "PRIORITY_4_NETWORK",
        "PRIORITY_5_EXPLORATORY"
      )
    )
  )
)

priority_order <- c(
  "PRIORITY_1_ROBUST",
  "PRIORITY_2_MULTI_LAYER",
  "PRIORITY_3_TRANSCRIPTOMIC",
  "PRIORITY_4_NETWORK",
  "PRIORITY_5_EXPLORATORY"
)

integrated$final_priority <- factor(
  integrated$final_priority,
  levels = priority_order
)

integrated <- integrated[
  order(
    integrated$final_priority,
    -integrated$validation_score,
    integrated$validation_rank
  ),
  ,
  drop = FALSE
]

rownames(integrated) <- NULL

integrated$final_priority_rank <- seq_len(nrow(integrated))


# ------------------------------------------------------------
# 8. ROBUST CANDIDATE PROFILE
# ------------------------------------------------------------

robust <- integrated[
  integrated$validation_class == "ROBUST_CROSS_LAYER",
  ,
  drop = FALSE
]

robust <- robust[
  order(-robust$validation_score),
  ,
  drop = FALSE
]

robust_profile <- robust[
  ,
  c(
    "final_priority_rank",
    "gene_symbol",
    "validation_score",
    "network_score",
    "logFC_siNRP1_3",
    "FDR_siNRP1_3",
    "logFC_siNRP1_5",
    "FDR_siNRP1_5",
    "direction_concordance",
    "expression_pattern",
    "pathway_count",
    "pathway_list",
    "module_count",
    "module_list"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 9. FINAL SHORTLIST PROFILE
# ------------------------------------------------------------

shortlist <- integrated[
  integrated$validation_class %in% c(
    "ROBUST_CROSS_LAYER",
    "MULTI_LAYER_SUPPORTED"
  ),
  ,
  drop = FALSE
]

shortlist_profile <- shortlist[
  ,
  c(
    "final_priority_rank",
    "gene_symbol",
    "validation_score",
    "validation_class",
    "network_score",
    "logFC_siNRP1_3",
    "FDR_siNRP1_3",
    "logFC_siNRP1_5",
    "FDR_siNRP1_5",
    "direction_concordance",
    "expression_pattern",
    "pathway_count",
    "pathway_list",
    "module_count",
    "module_list"
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 10. FOCUSED BIOLOGICAL PANEL
# ------------------------------------------------------------

focused_genes <- c(
  "THBS1",
  "CCN3",
  "TGFB2",
  "BMP4",
  "CAV1",
  "NRP1",
  "LYN",
  "TGFBR3",
  "ENG",
  "TGFBR2",
  "CCN2"
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

focused_panel$biological_interpretation <- ifelse(
  focused_panel$expression_pattern == "DOWN_BOTH",
  "Concordant decrease in both knockdown contrasts",
  ifelse(
    focused_panel$expression_pattern == "UP_BOTH",
    "Concordant increase in both knockdown contrasts",
    ifelse(
      focused_panel$expression_pattern == "OPPOSITE_DIRECTIONS",
      "Opposite expression directions across contrasts",
      ifelse(
        focused_panel$expression_pattern ==
          "SIGNIFICANT_ONLY_siNRP1_3",
        "Significant only in siNRP1.3",
        ifelse(
          focused_panel$expression_pattern ==
            "SIGNIFICANT_ONLY_siNRP1_5",
          "Significant only in siNRP1.5",
          "No significant differential expression in both contrasts"
        )
      )
    )
  )
)


# ------------------------------------------------------------
# 11. EXPRESSION PATTERN SUMMARY
# ------------------------------------------------------------

pattern_summary <- as.data.frame(
  table(
    factor(
      shortlist$expression_pattern,
      levels = c(
        "DOWN_BOTH",
        "UP_BOTH",
        "OPPOSITE_DIRECTIONS",
        "ZERO_LOGFC_IN_ONE_CONTRAST",
        "SIGNIFICANT_ONLY_siNRP1_3",
        "SIGNIFICANT_ONLY_siNRP1_5",
        "NOT_SIGNIFICANT"
      )
    )
  ),
  stringsAsFactors = FALSE
)

names(pattern_summary) <- c(
  "expression_pattern",
  "candidate_count"
)

pattern_summary$percentage <- if (nrow(shortlist) > 0) {
  100 * pattern_summary$candidate_count / nrow(shortlist)
} else {
  0
}


# ------------------------------------------------------------
# 12. PATHWAY FREQUENCY
# ------------------------------------------------------------

cat("Summarizing pathway evidence...\n")

pathway_frequency <- annotation_frequency(
  shortlist,
  "pathway_list"
)

if (nrow(pathway_frequency) > 0) {
  pathway_frequency <- pathway_frequency[
    order(-pathway_frequency$candidate_count),
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------
# 13. MODULE FREQUENCY
# ------------------------------------------------------------

cat("Summarizing module evidence...\n")

module_frequency <- annotation_frequency(
  shortlist,
  "module_list"
)

if (nrow(module_frequency) > 0) {
  module_frequency <- module_frequency[
    order(-module_frequency$candidate_count),
    ,
    drop = FALSE
  ]
}


# ------------------------------------------------------------
# 14. ROBUST CANDIDATE EXPRESSION SUMMARY
# ------------------------------------------------------------

robust_pattern_summary <- as.data.frame(
  table(
    factor(
      robust$expression_pattern,
      levels = c(
        "DOWN_BOTH",
        "UP_BOTH",
        "OPPOSITE_DIRECTIONS",
        "ZERO_LOGFC_IN_ONE_CONTRAST"
      )
    )
  ),
  stringsAsFactors = FALSE
)

names(robust_pattern_summary) <- c(
  "expression_pattern",
  "candidate_count"
)


# ------------------------------------------------------------
# 15. CANDIDATE EVIDENCE MATRIX
# ------------------------------------------------------------

evidence_matrix <- data.frame(
  gene_symbol = shortlist$gene_symbol,
  validation_class = shortlist$validation_class,
  final_priority = as.character(shortlist$final_priority),
  validation_score = shortlist$validation_score,
  significant_siNRP1_3 = shortlist$significant_siNRP1_3,
  significant_siNRP1_5 = shortlist$significant_siNRP1_5,
  significant_both = shortlist$significant_both,
  direction_concordance = shortlist$direction_concordance,
  expression_pattern = shortlist$expression_pattern,
  network_evidence = shortlist$has_network_evidence,
  pathway_evidence = shortlist$has_pathway_evidence,
  module_evidence = shortlist$has_module_evidence,
  pathway_count = shortlist$pathway_count,
  module_count = shortlist$module_count,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 16. CLASS DISTRIBUTION
# ------------------------------------------------------------

class_summary <- as.data.frame(
  table(
    factor(
      integrated$validation_class,
      levels = c(
        "ROBUST_CROSS_LAYER",
        "MULTI_LAYER_SUPPORTED",
        "TRANSCRIPTOMICALLY_SUPPORTED",
        "NETWORK_PRIORITY_ONLY",
        "EXPLORATORY"
      )
    )
  ),
  stringsAsFactors = FALSE
)

names(class_summary) <- c(
  "validation_class",
  "gene_count"
)

class_summary$percentage <- (
  100 * class_summary$gene_count / nrow(integrated)
)


# ------------------------------------------------------------
# 17. KEY CANDIDATE SUMMARY
# ------------------------------------------------------------

key_candidates <- c(
  "THBS1",
  "CCN3",
  "NRP1",
  "TGFBR3",
  "VEGFC",
  "PIK3CA",
  "FGF1",
  "RELN",
  "MET",
  "PAK2",
  "CCN2",
  "TGFB2",
  "BMP4",
  "CAV1",
  "LYN",
  "ENG",
  "TGFBR2"
)

key_candidate_table <- integrated[
  integrated$gene_symbol %in% key_candidates,
  ,
  drop = FALSE
]

key_candidate_table <- key_candidate_table[
  order(-key_candidate_table$validation_score),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 18. TOP PATHWAYS AND MODULES
# ------------------------------------------------------------

top_pathways <- head(pathway_frequency, 15)
top_modules <- head(module_frequency, 15)


# ------------------------------------------------------------
# 19. FINAL TABLE EXPORT
# ------------------------------------------------------------

cat("Exporting final candidate tables...\n")

write_csv(
  integrated,
  file.path(
    table_dir,
    "final_integrated_candidate_ranking.csv"
  )
)

write_csv(
  shortlist_profile,
  file.path(
    table_dir,
    "final_candidate_shortlist.csv"
  )
)

write_csv(
  robust_profile,
  file.path(
    table_dir,
    "robust_cross_layer_candidate_profile.csv"
  )
)

write_csv(
  focused_panel,
  file.path(
    table_dir,
    "focused_biological_candidate_panel.csv"
  )
)

write_csv(
  evidence_matrix,
  file.path(
    table_dir,
    "final_candidate_evidence_matrix.csv"
  )
)

write_csv(
  pattern_summary,
  file.path(
    table_dir,
    "candidate_expression_pattern_summary.csv"
  )
)

write_csv(
  robust_pattern_summary,
  file.path(
    table_dir,
    "robust_candidate_expression_summary.csv"
  )
)

write_csv(
  pathway_frequency,
  file.path(
    table_dir,
    "shortlist_pathway_frequency.csv"
  )
)

write_csv(
  module_frequency,
  file.path(
    table_dir,
    "shortlist_module_frequency.csv"
  )
)

write_csv(
  class_summary,
  file.path(
    table_dir,
    "final_candidate_class_distribution.csv"
  )
)

write_csv(
  key_candidate_table,
  file.path(
    table_dir,
    "key_biological_candidate_summary.csv"
  )
)


# ------------------------------------------------------------
# 20. BIOLOGICAL INTERPRETATION REPORT
# ------------------------------------------------------------

cat("Generating biological interpretation report...\n")

down_both <- shortlist$gene_symbol[
  shortlist$expression_pattern == "DOWN_BOTH"
]

up_both <- shortlist$gene_symbol[
  shortlist$expression_pattern == "UP_BOTH"
]

opposite <- shortlist$gene_symbol[
  shortlist$expression_pattern == "OPPOSITE_DIRECTIONS"
]

only_3_genes <- shortlist$gene_symbol[
  shortlist$expression_pattern ==
    "SIGNIFICANT_ONLY_siNRP1_3"
]

only_5_genes <- shortlist$gene_symbol[
  shortlist$expression_pattern ==
    "SIGNIFICANT_ONLY_siNRP1_5"
]

robust_down <- robust$gene_symbol[
  robust$expression_pattern == "DOWN_BOTH"
]

robust_up <- robust$gene_symbol[
  robust$expression_pattern == "UP_BOTH"
]

robust_opposite <- robust$gene_symbol[
  robust$expression_pattern == "OPPOSITE_DIRECTIONS"
]

top10 <- head(shortlist, 10)

top10_lines <- if (nrow(top10) > 0) {
  sprintf(
    "%d. %s | Score %.2f | %s | %s",
    seq_len(nrow(top10)),
    top10$gene_symbol,
    top10$validation_score,
    top10$validation_class,
    top10$expression_pattern
  )
} else {
  "No candidates available."
}

pathway_lines <- if (nrow(top_pathways) > 0) {
  sprintf(
    "%s | %d candidates",
    top_pathways$annotation,
    top_pathways$candidate_count
  )
} else {
  "No pathway annotations available."
}

module_lines <- if (nrow(top_modules) > 0) {
  sprintf(
    "%s | %d candidates",
    top_modules$annotation,
    top_modules$candidate_count
  )
} else {
  "No module annotations available."
}

report <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "SCRIPT 42 - FINAL BIOLOGICAL INTERPRETATION",
  "AND CANDIDATE PRIORITIZATION",
  "============================================================",
  "",
  "1. ANALYSIS OVERVIEW",
  "--------------------",
  paste("Integrated candidates:", nrow(integrated)),
  paste("Robust cross-layer candidates:", nrow(robust)),
  paste("Multi-layer supported candidates:", nrow(multi_layer)),
  paste("Final shortlist:", nrow(shortlist)),
  "",
  "The final prioritization is based on the integrated results",
  "generated by Script 41. No new differential expression,",
  "pathway enrichment or network inference was performed.",
  "",
  "2. CANDIDATE PRIORITIZATION",
  "---------------------------",
  paste(
    "Priority 1 - Robust cross-layer:",
    nrow(robust)
  ),
  paste(
    "Priority 2 - Multi-layer supported:",
    nrow(multi_layer)
  ),
  paste(
    "Priority 3 - Transcriptomically supported:",
    nrow(transcriptomic)
  ),
  paste(
    "Priority 4 - Network priority only:",
    nrow(network_only)
  ),
  paste(
    "Priority 5 - Exploratory:",
    nrow(exploratory)
  ),
  "",
  "3. TOP 10 FINAL CANDIDATES",
  "--------------------------",
  top10_lines,
  "",
  "4. EXPRESSION PATTERNS IN THE FINAL SHORTLIST",
  "----------------------------------------------",
  paste(
    "Concordant decrease in both contrasts:",
    length(down_both)
  ),
  paste(
    "Concordant increase in both contrasts:",
    length(up_both)
  ),
  paste(
    "Opposite directions:",
    length(opposite)
  ),
  paste(
    "Significant only in siNRP1.3:",
    length(only_3_genes)
  ),
  paste(
    "Significant only in siNRP1.5:",
    length(only_5_genes)
  ),
  "",
  paste(
    "Genes decreased in both:",
    format_gene_list(down_both)
  ),
  "",
  paste(
    "Genes increased in both:",
    format_gene_list(up_both)
  ),
  "",
  paste(
    "Genes with opposite directions:",
    format_gene_list(opposite)
  ),
  "",
  "5. ROBUST CROSS-LAYER CANDIDATES",
  "--------------------------------",
  paste(
    "Robust candidates decreased in both:",
    length(robust_down)
  ),
  paste(
    "Robust candidates increased in both:",
    length(robust_up)
  ),
  paste(
    "Robust candidates with opposite directions:",
    length(robust_opposite)
  ),
  "",
  paste(
    "Decreased in both:",
    format_gene_list(robust_down)
  ),
  "",
  paste(
    "Increased in both:",
    format_gene_list(robust_up)
  ),
  "",
  paste(
    "Opposite directions:",
    format_gene_list(robust_opposite)
  ),
  "",
  "6. MOST REPRESENTED PATHWAY ANNOTATIONS",
  "---------------------------------------",
  pathway_lines,
  "",
  "7. MOST REPRESENTED MODULE ANNOTATIONS",
  "--------------------------------------",
  module_lines,
  "",
  "8. FOCUSED BIOLOGICAL CANDIDATES",
  "--------------------------------",
  paste(
    "Focused genes represented in the integrated ranking:",
    nrow(focused_panel)
  ),
  paste(
    "Focused genes classified as robust:",
    sum(
      focused_panel$validation_class ==
        "ROBUST_CROSS_LAYER"
    )
  ),
  "",
  "9. SCIENTIFIC INTERPRETATION",
  "----------------------------",
  "The integrated prioritization identifies a subset of genes",
  "supported by multiple evidence layers, including network",
  "ranking, differential expression, pathway annotations and",
  "gene-module assignments.",
  "",
  "THBS1, CCN3 and NRP1 occupy the leading positions in the",
  "previously generated integrated ranking. Their prioritization",
  "reflects the combination of transcriptomic and network-based",
  "evidence rather than a direct measurement of biological effect.",
  "",
  "The comparison between siNRP1.3 and siNRP1.5 distinguishes",
  "genes with concordant expression patterns from genes showing",
  "contrast-specific or opposite responses.",
  "",
  "Concordant expression changes across both knockdown conditions",
  "provide evidence of consistency across perturbation conditions.",
  "They do not independently establish a causal relationship.",
  "",
  "Pathway and module annotations provide functional context,",
  "but they are derived from overlapping biological knowledge",
  "and should not be interpreted as independent experiments.",
  "",
  "10. LIMITATIONS",
  "---------------",
  "The candidate score is a descriptive prioritization metric.",
  "It is not a probability, effect size or statistical test.",
  "",
  "The two knockdown conditions may differ in perturbation strength",
  "and should not be assumed to represent identical biological",
  "responses.",
  "",
  "Pathway and module support depends on annotation coverage.",
  "Absence of an annotation does not demonstrate absence of",
  "biological involvement.",
  "",
  "The results are based on bulk RNA sequencing and computational",
  "integration. Experimental validation is required to establish",
  "the functional relevance of prioritized candidates.",
  "",
  "11. CONCLUSION",
  "---------------",
  "The final candidate ranking provides a structured shortlist",
  "for downstream biological interpretation and experimental",
  "prioritization.",
  "",
  "The 23 robust cross-layer candidates represent the highest",
  "evidence-convergence group under the criteria defined in",
  "Script 41. The remaining 56 multi-layer supported candidates",
  "provide a secondary group for further investigation.",
  "",
  "============================================================",
  "END OF SCRIPT 42 REPORT",
  "============================================================"
)

writeLines(
  report,
  file.path(
    report_dir,
    "final_biological_interpretation_report.txt"
  )
)


# ------------------------------------------------------------
# 21. FINAL AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "integrated_candidates",
    "robust_cross_layer_candidates",
    "multi_layer_supported_candidates",
    "transcriptomically_supported_candidates",
    "network_priority_only_candidates",
    "exploratory_candidates",
    "final_shortlist",
    "shortlist_down_both",
    "shortlist_up_both",
    "shortlist_opposite_directions",
    "shortlist_significant_only_siNRP1_3",
    "shortlist_significant_only_siNRP1_5",
    "robust_down_both",
    "robust_up_both",
    "robust_opposite_directions",
    "focused_panel_genes",
    "shortlist_pathway_annotations",
    "shortlist_module_annotations"
  ),
  value = c(
    nrow(integrated),
    nrow(robust),
    nrow(multi_layer),
    nrow(transcriptomic),
    nrow(network_only),
    nrow(exploratory),
    nrow(shortlist),
    length(down_both),
    length(up_both),
    length(opposite),
    length(only_3_genes),
    length(only_5_genes),
    length(robust_down),
    length(robust_up),
    length(robust_opposite),
    nrow(focused_panel),
    nrow(pathway_frequency),
    nrow(module_frequency)
  ),
  stringsAsFactors = FALSE
)

write_csv(
  audit,
  file.path(audit_dir, "script42_audit.csv")
)


# ------------------------------------------------------------
# 22. FINAL VALIDATION
# ------------------------------------------------------------

if (sum(class_summary$gene_count) != nrow(integrated)) {
  stop("Class distribution does not match integrated candidates.")
}

if (nrow(shortlist) != nrow(robust) + nrow(multi_layer)) {
  stop("Shortlist size does not match candidate classes.")
}

if (anyDuplicated(shortlist$gene_symbol)) {
  stop("Duplicated genes in final shortlist.")
}

if (anyDuplicated(robust$gene_symbol)) {
  stop("Duplicated genes in robust candidate group.")
}

if (any(
  robust$validation_class != "ROBUST_CROSS_LAYER"
)) {
  stop("Invalid candidate in robust group.")
}

if (any(
  shortlist$validation_class %in% c(
    "TRANSCRIPTOMICALLY_SUPPORTED",
    "NETWORK_PRIORITY_ONLY",
    "EXPLORATORY"
  )
)) {
  stop("Unexpected candidate class in final shortlist.")
}

if (sum(pattern_summary$candidate_count) != nrow(shortlist)) {
  stop("Expression pattern counts do not match shortlist.")
}

# Verify that robust candidates are significant in both contrasts.

if (any(
  !robust$significant_siNRP1_3 |
  !robust$significant_siNRP1_5
)) {
  stop("Robust candidate significance validation failed.")
}

# Verify that robust candidates have pathway and module support.

if (any(
  robust$pathway_count == 0 |
  robust$module_count == 0
)) {
  stop("Robust candidate annotation validation failed.")
}


# ------------------------------------------------------------
# 23. COMPLETION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 42 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("INTEGRATED CANDIDATES:", nrow(integrated), "\n")
cat("ROBUST CROSS-LAYER:", nrow(robust), "\n")
cat("MULTI-LAYER SUPPORTED:", nrow(multi_layer), "\n")
cat("FINAL SHORTLIST:", nrow(shortlist), "\n\n")

cat("EXPRESSION PATTERNS - FINAL SHORTLIST\n")
print(pattern_summary, row.names = FALSE)

cat("\nTOP 10 FINAL CANDIDATES\n")

print(
  top10[
    ,
    c(
      "gene_symbol",
      "validation_score",
      "validation_class",
      "expression_pattern"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nROBUST CROSS-LAYER CANDIDATES\n")

print(
  robust_profile[
    ,
    c(
      "gene_symbol",
      "validation_score",
      "expression_pattern",
      "pathway_count",
      "module_count"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nOUTPUT DIRECTORY:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
