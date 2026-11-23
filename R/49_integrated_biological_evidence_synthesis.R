
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 49 - INTEGRATED BIOLOGICAL EVIDENCE SYNTHESIS
#
# Objective:
# Build a traceable evidence framework for the 79 final
# candidate genes by integrating:
#   1. Differential expression
#   2. Network evidence
#   3. Pathway evidence
#   4. Functional module evidence
#
# Existing validation scores and candidate rankings are
# preserved. No new composite score is calculated.
#
# Input:
#   Script 37 - Gene-pathway and gene-module edges
#   Script 41 - Integrated candidate validation
#   Script 42 - Final biological interpretation
#
# Output:
#   results/integrated_biological_evidence_synthesis/
#
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 49 - INTEGRATED BIOLOGICAL EVIDENCE SYNTHESIS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

input_dir <- file.path(project_dir, "results")

evidence_dir <- file.path(
  input_dir,
  "final_biological_interpretation",
  "tables"
)

validation_dir <- file.path(
  input_dir,
  "integrated_candidate_gene_validation",
  "tables"
)

network_dir <- file.path(
  input_dir,
  "gene_pathway_module_network",
  "tables"
)

output_dir <- file.path(
  input_dir,
  "integrated_biological_evidence_synthesis"
)

tables_dir <- file.path(output_dir, "tables")
reports_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (directory in c(
  tables_dir,
  reports_dir,
  audit_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# 2. INPUT FILES
# ------------------------------------------------------------

input_files <- c(
  final_evidence_matrix = file.path(
    evidence_dir,
    "final_candidate_evidence_matrix.csv"
  ),
  
  final_candidate_ranking = file.path(
    evidence_dir,
    "final_integrated_candidate_ranking.csv"
  ),
  
  final_candidate_shortlist = file.path(
    evidence_dir,
    "final_candidate_shortlist.csv"
  ),
  
  integrated_validation = file.path(
    validation_dir,
    "integrated_candidate_gene_validation.csv"
  ),
  
  gene_pathway_edges = file.path(
    network_dir,
    "gene_pathway_edges.csv"
  ),
  
  gene_module_edges = file.path(
    network_dir,
    "gene_module_edges.csv"
  )
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required input files:",
      paste(missing_files, collapse = "\n"),
      sep = "\n"
    )
  )
}

# ------------------------------------------------------------
# 3. LOGGING
# ------------------------------------------------------------

log_file <- file.path(
  audit_dir,
  "script49.log"
)

log_connection <- file(
  log_file,
  open = "wt"
)

sink(log_connection, split = TRUE)

cat("Project:", project_dir, "\n")
cat("Execution date:", as.character(Sys.time()), "\n\n")

# ------------------------------------------------------------
# 4. HELPER FUNCTIONS
# ------------------------------------------------------------

read_input <- function(path) {
  
  data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  
  names(data) <- trimws(names(data))
  
  data
}

check_columns <- function(data, required, label) {
  
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing columns in ",
        label,
        ": ",
        paste(missing, collapse = ", ")
      )
    )
  }
}

to_boolean <- function(x) {
  
  value <- toupper(trimws(as.character(x)))
  
  result <- rep(NA, length(value))
  
  result[value %in% c(
    "TRUE", "T", "1", "YES", "Y"
  )] <- TRUE
  
  result[value %in% c(
    "FALSE", "F", "0", "NO", "N"
  )] <- FALSE
  
  result
}

to_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(as.character(x))
  )
}

collapse_unique <- function(x, separator = " | ") {
  
  values <- unique(as.character(x))
  
  values <- values[
    !is.na(values) &
      nzchar(trimws(values))
  ]
  
  if (length(values) == 0) {
    return(NA_character_)
  }
  
  paste(sort(values), collapse = separator)
}

write_table <- function(data, filename) {
  
  write.csv(
    data,
    file.path(tables_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

summarise_edges <- function(
    data,
    grouping_column,
    metadata_columns,
    gene_column = "gene_symbol") {
  
  if (nrow(data) == 0) {
    return(data.frame())
  }
  
  groups <- split(
    data,
    as.character(data[[grouping_column]])
  )
  
  results <- lapply(groups, function(group_data) {
    
    first_row <- group_data[1, metadata_columns, drop = FALSE]
    
    first_row$candidate_gene_count <- length(
      unique(group_data[[gene_column]])
    )
    
    first_row$candidate_genes <- collapse_unique(
      group_data[[gene_column]],
      separator = ";"
    )
    
    first_row
  })
  
  result <- do.call(rbind, results)
  
  rownames(result) <- NULL
  
  result <- result[
    order(
      -result$candidate_gene_count,
      as.character(result[[grouping_column]])
    ),
    ,
    drop = FALSE
  ]
  
  result
}

# ------------------------------------------------------------
# 5. LOAD DATA
# ------------------------------------------------------------

cat("STEP 1 - Loading input datasets\n")

evidence_matrix <- read_input(
  input_files["final_evidence_matrix"]
)

candidate_ranking <- read_input(
  input_files["final_candidate_ranking"]
)

candidate_shortlist <- read_input(
  input_files["final_candidate_shortlist"]
)

integrated_validation <- read_input(
  input_files["integrated_validation"]
)

pathway_edges <- read_input(
  input_files["gene_pathway_edges"]
)

module_edges <- read_input(
  input_files["gene_module_edges"]
)

cat("Evidence matrix:", nrow(evidence_matrix), "rows\n")
cat("Candidate ranking:", nrow(candidate_ranking), "rows\n")
cat("Final shortlist:", nrow(candidate_shortlist), "rows\n")
cat("Integrated validation:", nrow(integrated_validation), "rows\n")
cat("Gene-pathway edges:", nrow(pathway_edges), "rows\n")
cat("Gene-module edges:", nrow(module_edges), "rows\n\n")

# ------------------------------------------------------------
# 6. INPUT VALIDATION
# ------------------------------------------------------------

cat("STEP 2 - Validating input structure\n")

check_columns(
  evidence_matrix,
  c(
    "gene_symbol",
    "validation_class",
    "final_priority",
    "validation_score"
  ),
  "final evidence matrix"
)

check_columns(
  candidate_ranking,
  c(
    "gene_symbol",
    "validation_score",
    "validation_class",
    "final_priority",
    "logFC_siNRP1_3",
    "FDR_siNRP1_3",
    "logFC_siNRP1_5",
    "FDR_siNRP1_5",
    "has_network_evidence",
    "has_DE_evidence",
    "has_pathway_evidence",
    "has_module_evidence",
    "evidence_layer_count"
  ),
  "final candidate ranking"
)

check_columns(
  candidate_shortlist,
  c(
    "gene_symbol",
    "final_priority_rank",
    "validation_score"
  ),
  "final candidate shortlist"
)

check_columns(
  integrated_validation,
  c(
    "gene_symbol",
    "validation_score",
    "validation_class"
  ),
  "integrated validation"
)

check_columns(
  pathway_edges,
  c(
    "pathway_id",
    "gene_symbol",
    "term",
    "ontology",
    "group"
  ),
  "gene-pathway edges"
)

check_columns(
  module_edges,
  c(
    "gene_symbol",
    "module",
    "pathway_id"
  ),
  "gene-module edges"
)

for (dataset_name in c(
  "evidence_matrix",
  "candidate_ranking",
  "candidate_shortlist",
  "integrated_validation",
  "pathway_edges",
  "module_edges"
)) {
  
  dataset <- get(dataset_name)
  
  if (any(is.na(dataset$gene_symbol)) ||
      any(!nzchar(trimws(dataset$gene_symbol)))) {
    
    stop(
      paste(
        "Invalid gene symbols detected in",
        dataset_name
      )
    )
  }
}

if (anyDuplicated(candidate_shortlist$gene_symbol)) {
  stop("Duplicate genes detected in final candidate shortlist.")
}

if (anyDuplicated(candidate_ranking$gene_symbol)) {
  stop("Duplicate genes detected in final candidate ranking.")
}

if (anyDuplicated(evidence_matrix$gene_symbol)) {
  stop("Duplicate genes detected in final evidence matrix.")
}

candidate_genes <- unique(
  candidate_shortlist$gene_symbol
)

cat(
  "Unique final candidates:",
  length(candidate_genes),
  "\n"
)

# ------------------------------------------------------------
# 7. BUILD INTEGRATED CANDIDATE EVIDENCE LEDGER
# ------------------------------------------------------------

cat("\nSTEP 3 - Building integrated evidence ledger\n")

# Preserve the original shortlist and its ranking.

ledger <- candidate_shortlist

# Add evidence descriptions from Script 42.

evidence_columns <- intersect(
  c(
    "gene_symbol",
    "network_evidence",
    "pathway_evidence",
    "module_evidence",
    "pathway_count",
    "module_count"
  ),
  names(evidence_matrix)
)

evidence_additional <- evidence_matrix[
  ,
  evidence_columns,
  drop = FALSE
]

evidence_additional <- evidence_additional[
  !duplicated(evidence_additional$gene_symbol),
  ,
  drop = FALSE
]

ledger <- merge(
  ledger,
  evidence_additional,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_evidence")
)

# Add evidence-layer indicators and component scores.

ranking_columns <- intersect(
  c(
    "gene_symbol",
    "network_class",
    "evidence_layer_count",
    "has_network_evidence",
    "has_DE_evidence",
    "has_pathway_evidence",
    "has_module_evidence",
    "network_component",
    "DE_component",
    "concordance_component",
    "pathway_component",
    "module_component",
    "direction_siNRP1_3",
    "direction_siNRP1_5",
    "significant_both",
    "significant_either",
    "significant_siNRP1_3",
    "significant_siNRP1_5"
  ),
  names(candidate_ranking)
)

ranking_additional <- candidate_ranking[
  ,
  ranking_columns,
  drop = FALSE
]

ranking_additional <- ranking_additional[
  !duplicated(ranking_additional$gene_symbol),
  ,
  drop = FALSE
]

ledger <- merge(
  ledger,
  ranking_additional,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

# Preserve the original candidate order.

ledger$original_order <- match(
  ledger$gene_symbol,
  candidate_genes
)

ledger <- ledger[
  order(ledger$original_order),
  ,
  drop = FALSE
]

ledger$original_order <- NULL

# Convert evidence indicators to logical values.

evidence_flags <- c(
  "has_network_evidence",
  "has_DE_evidence",
  "has_pathway_evidence",
  "has_module_evidence",
  "significant_both",
  "significant_either",
  "significant_siNRP1_3",
  "significant_siNRP1_5"
)

for (column in intersect(evidence_flags, names(ledger))) {
  ledger[[column]] <- to_boolean(ledger[[column]])
}

cat(
  "Integrated candidate ledger:",
  nrow(ledger),
  "genes\n"
)

# ------------------------------------------------------------
# 8. CANDIDATE PATHWAY EVIDENCE
# ------------------------------------------------------------

cat("\nSTEP 4 - Extracting candidate pathway evidence\n")

candidate_pathway_edges <- pathway_edges[
  pathway_edges$gene_symbol %in% candidate_genes,
  ,
  drop = FALSE
]

candidate_pathway_edges <- unique(
  candidate_pathway_edges
)

candidate_pathway_edges <- candidate_pathway_edges[
  order(
    candidate_pathway_edges$gene_symbol,
    candidate_pathway_edges$ontology,
    candidate_pathway_edges$term
  ),
  ,
  drop = FALSE
]

cat(
  "Candidate pathway associations:",
  nrow(candidate_pathway_edges),
  "\n"
)

# Pathway support among final candidates.

pathway_summary <- summarise_edges(
  candidate_pathway_edges,
  grouping_column = "pathway_id",
  metadata_columns = c(
    "pathway_id",
    "term",
    "ontology",
    "group"
  )
)

if (nrow(pathway_summary) > 0) {
  
  pathway_summary <- pathway_summary[
    order(
      -pathway_summary$candidate_gene_count,
      pathway_summary$term
    ),
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------
# 9. CANDIDATE MODULE EVIDENCE
# ------------------------------------------------------------

cat("\nSTEP 5 - Extracting candidate module evidence\n")

candidate_module_edges <- module_edges[
  module_edges$gene_symbol %in% candidate_genes,
  ,
  drop = FALSE
]

candidate_module_edges <- unique(
  candidate_module_edges
)

candidate_module_edges <- candidate_module_edges[
  order(
    candidate_module_edges$module,
    candidate_module_edges$gene_symbol
  ),
  ,
  drop = FALSE
]

cat(
  "Candidate module associations:",
  nrow(candidate_module_edges),
  "\n"
)

module_summary <- summarise_edges(
  candidate_module_edges,
  grouping_column = "module",
  metadata_columns = c("module")
)

if (nrow(module_summary) > 0) {
  
  module_summary <- module_summary[
    order(
      -module_summary$candidate_gene_count,
      module_summary$module
    ),
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------
# 10. CANDIDATE-LEVEL PATHWAY AND MODULE COUNTS
# ------------------------------------------------------------

cat("\nSTEP 6 - Calculating candidate evidence coverage\n")

if (nrow(candidate_pathway_edges) > 0) {
  
  pathway_counts <- aggregate(
    pathway_id ~ gene_symbol,
    data = candidate_pathway_edges,
    FUN = function(x) length(unique(x))
  )
  
  names(pathway_counts)[2] <- "observed_pathway_count"
  
} else {
  
  pathway_counts <- data.frame(
    gene_symbol = character(),
    observed_pathway_count = integer()
  )
}

if (nrow(candidate_module_edges) > 0) {
  
  module_counts <- aggregate(
    module ~ gene_symbol,
    data = candidate_module_edges,
    FUN = function(x) length(unique(x))
  )
  
  names(module_counts)[2] <- "observed_module_count"
  
} else {
  
  module_counts <- data.frame(
    gene_symbol = character(),
    observed_module_count = integer()
  )
}

ledger <- merge(
  ledger,
  pathway_counts,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

ledger <- merge(
  ledger,
  module_counts,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

ledger$observed_pathway_count[
  is.na(ledger$observed_pathway_count)
] <- 0

ledger$observed_module_count[
  is.na(ledger$observed_module_count)
] <- 0

ledger$observed_pathway_count <- as.integer(
  ledger$observed_pathway_count
)

ledger$observed_module_count <- as.integer(
  ledger$observed_module_count
)

ledger$original_order <- match(
  ledger$gene_symbol,
  candidate_genes
)

ledger <- ledger[
  order(ledger$original_order),
  ,
  drop = FALSE
]

ledger$original_order <- NULL

# ------------------------------------------------------------
# 11. EVIDENCE LAYER COVERAGE
# ------------------------------------------------------------

cat("\nSTEP 7 - Summarizing evidence layer coverage\n")

layer_columns <- c(
  "has_DE_evidence",
  "has_network_evidence",
  "has_pathway_evidence",
  "has_module_evidence"
)

available_layer_columns <- intersect(
  layer_columns,
  names(ledger)
)

if (length(available_layer_columns) != 4) {
  stop(
    "Not all four evidence-layer indicators are available."
  )
}

layer_matrix <- as.data.frame(
  lapply(
    ledger[, available_layer_columns, drop = FALSE],
    function(x) {
      ifelse(is.na(x), FALSE, as.logical(x))
    }
  )
)

names(layer_matrix) <- c(
  "DE",
  "Network",
  "Pathway",
  "Module"
)

ledger$observed_layer_count <- rowSums(
  layer_matrix
)

layer_coverage <- data.frame(
  evidence_layer = c(
    "Differential expression",
    "Network",
    "Pathway",
    "Functional module"
  ),
  candidate_count = c(
    sum(layer_matrix$DE),
    sum(layer_matrix$Network),
    sum(layer_matrix$Pathway),
    sum(layer_matrix$Module)
  ),
  candidate_percentage = round(
    c(
      mean(layer_matrix$DE),
      mean(layer_matrix$Network),
      mean(layer_matrix$Pathway),
      mean(layer_matrix$Module)
    ) * 100,
    2
  )
)

convergence_summary <- as.data.frame(
  table(
    factor(
      ledger$observed_layer_count,
      levels = 0:4
    )
  )
)

names(convergence_summary) <- c(
  "evidence_layers",
  "candidate_count"
)

convergence_summary$evidence_layers <- as.integer(
  as.character(convergence_summary$evidence_layers)
)

convergence_summary$candidate_percentage <- round(
  convergence_summary$candidate_count /
    nrow(ledger) * 100,
  2
)

# ------------------------------------------------------------
# 12. DIRECTIONAL EXPRESSION CONSISTENCY
# ------------------------------------------------------------

cat("\nSTEP 8 - Summarizing expression direction\n")

if (all(c(
  "logFC_siNRP1_3",
  "logFC_siNRP1_5"
) %in% names(ledger))) {
  
  ledger$logFC_siNRP1_3 <- to_numeric(
    ledger$logFC_siNRP1_3
  )
  
  ledger$logFC_siNRP1_5 <- to_numeric(
    ledger$logFC_siNRP1_5
  )
  
  ledger$expression_direction_pair <- ifelse(
    is.na(ledger$logFC_siNRP1_3) |
      is.na(ledger$logFC_siNRP1_5),
    "NOT_AVAILABLE",
    ifelse(
      ledger$logFC_siNRP1_3 < 0 &
        ledger$logFC_siNRP1_5 < 0,
      "DOWNREGULATED_BOTH",
      ifelse(
        ledger$logFC_siNRP1_3 > 0 &
          ledger$logFC_siNRP1_5 > 0,
        "UPREGULATED_BOTH",
        ifelse(
          ledger$logFC_siNRP1_3 == 0 |
            ledger$logFC_siNRP1_5 == 0,
          "ZERO_LOGFC_IN_AT_LEAST_ONE_CONDITION",
          "OPPOSITE_DIRECTIONS"
        )
      )
    )
  )
  
  direction_summary <- as.data.frame(
    table(
      ledger$expression_direction_pair
    )
  )
  
  names(direction_summary) <- c(
    "expression_pattern",
    "candidate_count"
  )
  
  direction_summary$candidate_percentage <- round(
    direction_summary$candidate_count /
      nrow(ledger) * 100,
    2
  )
  
} else {
  
  direction_summary <- data.frame(
    expression_pattern = character(),
    candidate_count = integer(),
    candidate_percentage = numeric()
  )
}

# ------------------------------------------------------------
# 13. FINAL EVIDENCE TABLES
# ------------------------------------------------------------

cat("\nSTEP 9 - Exporting integrated evidence tables\n")

# Main candidate-level evidence ledger.

write_table(
  ledger,
  "integrated_candidate_evidence_ledger.csv"
)

# Candidate-to-pathway associations.

write_table(
  candidate_pathway_edges,
  "candidate_pathway_evidence.csv"
)

# Candidate-to-module associations.

write_table(
  candidate_module_edges,
  "candidate_module_evidence.csv"
)

# Pathway support among the final candidate panel.

write_table(
  pathway_summary,
  "candidate_pathway_support_summary.csv"
)

# Module support among the final candidate panel.

write_table(
  module_summary,
  "candidate_module_support_summary.csv"
)

# Evidence layer coverage.

write_table(
  layer_coverage,
  "candidate_evidence_layer_coverage.csv"
)

# Evidence convergence.

write_table(
  convergence_summary,
  "candidate_evidence_convergence_summary.csv"
)

# Expression direction consistency.

write_table(
  direction_summary,
  "candidate_expression_direction_summary.csv"
)

# ------------------------------------------------------------
# 14. PRIORITY-ORDERED CANDIDATE PANEL
# ------------------------------------------------------------

priority_panel <- ledger

if ("final_priority_rank" %in% names(priority_panel)) {
  
  priority_panel$final_priority_rank <- to_numeric(
    priority_panel$final_priority_rank
  )
  
  priority_panel <- priority_panel[
    order(
      priority_panel$final_priority_rank,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}

write_table(
  priority_panel,
  "priority_ordered_integrated_candidate_panel.csv"
)

# ------------------------------------------------------------
# 15. AUDIT
# ------------------------------------------------------------

cat("\nSTEP 10 - Generating audit\n")

audit <- data.frame(
  metric = c(
    "input_evidence_matrix_rows",
    "input_ranking_rows",
    "input_shortlist_rows",
    "input_integrated_validation_rows",
    "input_pathway_edge_rows",
    "input_module_edge_rows",
    "unique_final_candidates",
    "integrated_ledger_rows",
    "candidate_pathway_associations",
    "candidate_module_associations",
    "pathways_supported_by_candidates",
    "modules_supported_by_candidates",
    "candidates_with_DE_evidence",
    "candidates_with_network_evidence",
    "candidates_with_pathway_evidence",
    "candidates_with_module_evidence",
    "candidates_with_all_four_layers"
  ),
  value = c(
    nrow(evidence_matrix),
    nrow(candidate_ranking),
    nrow(candidate_shortlist),
    nrow(integrated_validation),
    nrow(pathway_edges),
    nrow(module_edges),
    length(candidate_genes),
    nrow(ledger),
    nrow(candidate_pathway_edges),
    nrow(candidate_module_edges),
    nrow(pathway_summary),
    nrow(module_summary),
    sum(layer_matrix$DE),
    sum(layer_matrix$Network),
    sum(layer_matrix$Pathway),
    sum(layer_matrix$Module),
    sum(ledger$observed_layer_count == 4)
  )
)

write.csv(
  audit,
  file.path(
    audit_dir,
    "script49_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. SCIENTIFIC REPORT
# ------------------------------------------------------------

cat("\nSTEP 11 - Generating synthesis report\n")

report_file <- file.path(
  reports_dir,
  "integrated_biological_evidence_synthesis_report.txt"
)

report_lines <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "SCRIPT 49 - INTEGRATED BIOLOGICAL EVIDENCE SYNTHESIS",
  "============================================================",
  "",
  paste("Execution date:", Sys.time()),
  "",
  "1. INPUT DATASETS",
  paste("Final evidence matrix:", nrow(evidence_matrix), "rows"),
  paste("Integrated candidate ranking:", nrow(candidate_ranking), "rows"),
  paste("Final candidate shortlist:", nrow(candidate_shortlist), "rows"),
  paste("Integrated validation:", nrow(integrated_validation), "rows"),
  paste("Gene-pathway associations:", nrow(pathway_edges), "rows"),
  paste("Gene-module associations:", nrow(module_edges), "rows"),
  "",
  "2. FINAL CANDIDATE PANEL",
  paste("Unique candidates:", length(candidate_genes)),
  paste("Evidence ledger records:", nrow(ledger)),
  "",
  "3. EVIDENCE LAYER COVERAGE",
  paste(
    layer_coverage$evidence_layer,
    ":",
    layer_coverage$candidate_count,
    "candidates (",
    layer_coverage$candidate_percentage,
    "%)",
    sep = ""
  ),
  "",
  "4. EVIDENCE CONVERGENCE"
)

for (i in seq_len(nrow(convergence_summary))) {
  
  report_lines <- c(
    report_lines,
    paste(
      convergence_summary$evidence_layers[i],
      "evidence layers:",
      convergence_summary$candidate_count[i],
      "candidates (",
      convergence_summary$candidate_percentage[i],
      "%)"
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "5. PATHWAY SUPPORT",
  paste(
    "Distinct pathways represented among final candidates:",
    nrow(pathway_summary)
  ),
  "",
  "Top pathways by number of represented candidates:"
)

if (nrow(pathway_summary) > 0) {
  
  top_pathways <- head(pathway_summary, 15)
  
  for (i in seq_len(nrow(top_pathways))) {
    
    report_lines <- c(
      report_lines,
      paste(
        top_pathways$term[i],
        "|",
        top_pathways$ontology[i],
        "| candidates:",
        top_pathways$candidate_gene_count[i],
        "| genes:",
        top_pathways$candidate_genes[i]
      )
    )
  }
}

report_lines <- c(
  report_lines,
  "",
  "6. FUNCTIONAL MODULE SUPPORT",
  paste(
    "Distinct modules represented among final candidates:",
    nrow(module_summary)
  ),
  "",
  "Module support:"
)

if (nrow(module_summary) > 0) {
  
  for (i in seq_len(nrow(module_summary))) {
    
    report_lines <- c(
      report_lines,
      paste(
        module_summary$module[i],
        "| candidates:",
        module_summary$candidate_gene_count[i],
        "| genes:",
        module_summary$candidate_genes[i]
      )
    )
  }
}

report_lines <- c(
  report_lines,
  "",
  "7. EXPRESSION DIRECTION",
  paste(
    direction_summary$expression_pattern,
    ":",
    direction_summary$candidate_count,
    "candidates"
  ),
  "",
  "8. TOP PRIORITY CANDIDATES"
)

top_candidates <- head(priority_panel, 15)

for (i in seq_len(nrow(top_candidates))) {
  
  report_lines <- c(
    report_lines,
    paste(
      top_candidates$final_priority_rank[i],
      top_candidates$gene_symbol[i],
      "| score:",
      round(
        to_numeric(top_candidates$validation_score[i]),
        3
      ),
      "| class:",
      top_candidates$validation_class[i]
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "9. INTERPRETATION AND LIMITATIONS",
  "",
  paste(
    "This synthesis integrates previously generated",
    "transcriptomic, network, pathway and module evidence."
  ),
  paste(
    "Existing validation scores and candidate rankings",
    "were preserved and were not recalculated."
  ),
  paste(
    "Pathway and module association counts represent",
    "the annotations present in the input network tables."
  ),
  paste(
    "They do not independently establish pathway activation,",
    "causal relationships or direct molecular interactions."
  ),
  paste(
    "Evidence convergence indicates support across analytical",
    "layers and should not be interpreted as independent",
    "experimental validation."
  ),
  paste(
    "The two siNRP1 conditions represent concordance within",
    "the analyzed dataset, not independent external replication."
  ),
  "",
  "============================================================",
  "END OF REPORT",
  "============================================================"
)

writeLines(
  report_lines,
  report_file
)

# ------------------------------------------------------------
# 17. FINAL CHECKS
# ------------------------------------------------------------

cat("\nSTEP 12 - Final consistency checks\n")

checks <- data.frame(
  check = c(
    "Shortlist and ledger have identical candidate counts",
    "No duplicated genes in evidence ledger",
    "All shortlisted genes represented in ledger",
    "All pathway associations belong to shortlisted genes",
    "All module associations belong to shortlisted genes",
    "Evidence layer count is between zero and four"
  ),
  passed = c(
    nrow(ledger) == length(candidate_genes),
    
    !anyDuplicated(ledger$gene_symbol),
    
    setequal(
      ledger$gene_symbol,
      candidate_genes
    ),
    
    all(
      candidate_pathway_edges$gene_symbol %in%
        candidate_genes
    ),
    
    all(
      candidate_module_edges$gene_symbol %in%
        candidate_genes
    ),
    
    all(
      ledger$observed_layer_count >= 0 &
        ledger$observed_layer_count <= 4
    )
  )
)

write.csv(
  checks,
  file.path(
    audit_dir,
    "script49_consistency_checks.csv"
  ),
  row.names = FALSE
)

print(checks)

if (!all(checks$passed)) {
  stop("One or more final consistency checks failed.")
}

# ------------------------------------------------------------
# 18. COMPLETION SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SCRIPT 49 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("Final candidates:", nrow(ledger), "\n")
cat("Pathway associations:", nrow(candidate_pathway_edges), "\n")
cat("Module associations:", nrow(candidate_module_edges), "\n")
cat("Distinct candidate-supported pathways:", nrow(pathway_summary), "\n")
cat("Distinct candidate-supported modules:", nrow(module_summary), "\n")
cat(
  "Candidates supported by all four evidence layers:",
  sum(ledger$observed_layer_count == 4),
  "\n"
)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nReport:\n")
cat(report_file, "\n")

cat("\nAudit:\n")
cat(file.path(audit_dir, "script49_audit.csv"), "\n")

cat("\n============================================================\n")

sink()
close(log_connection)
