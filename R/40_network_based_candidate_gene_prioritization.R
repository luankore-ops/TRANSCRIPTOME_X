
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 40 - NETWORK-BASED CANDIDATE GENE PRIORITIZATION
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 40 - NETWORK-BASED CANDIDATE GENE PRIORITIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = TRUE
)

network_dir <- file.path(
  project_dir, "results", "gene_pathway_module_network"
)

topology_dir <- file.path(
  project_dir, "results", "network_topology_analysis"
)

output_dir <- file.path(
  project_dir, "results", "network_candidate_prioritization"
)

table_dir <- file.path(output_dir, "tables")
report_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (d in c(table_dir, report_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# 2. INPUT FILES
# ------------------------------------------------------------

required_files <- c(
  gene_pathway_edges = file.path(
    network_dir, "tables", "gene_pathway_edges.csv"
  ),
  gene_module_edges = file.path(
    network_dir, "tables", "gene_module_edges.csv"
  ),
  node_topology = file.path(
    topology_dir, "tables", "node_topology_metrics.csv"
  )
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------

read_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

find_column <- function(data, candidates,
                        required = FALSE,
                        label = "data") {
  
  normalized <- tolower(
    gsub("[^a-z0-9]", "", names(data))
  )
  
  candidates_normalized <- tolower(
    gsub("[^a-z0-9]", "", candidates)
  )
  
  match_index <- match(
    candidates_normalized,
    normalized
  )
  
  match_index <- match_index[!is.na(match_index)]
  
  if (length(match_index) == 0) {
    if (required) {
      stop(
        "Cannot identify required column in ", label,
        ". Available columns: ",
        paste(names(data), collapse = ", ")
      )
    }
    return(NULL)
  }
  
  names(data)[match_index[1]]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

scale_rank <- function(x) {
  
  x <- safe_numeric(x)
  result <- rep(0, length(x))
  valid <- is.finite(x)
  
  if (sum(valid) == 0) return(result)
  
  if (sum(valid) == 1) {
    result[valid] <- 1
    return(result)
  }
  
  ranks <- rank(x[valid], ties.method = "average")
  
  result[valid] <- (
    ranks - 1
  ) / (
    length(ranks) - 1
  )
  
  result
}

count_unique <- function(data, column) {
  
  values <- clean_gene(data[[column]])
  values <- values[!is.na(values)]
  
  if (length(values) == 0) {
    return(data.frame(
      gene_symbol = character(),
      count = integer()
    ))
  }
  
  result <- as.data.frame(
    table(values),
    stringsAsFactors = FALSE
  )
  
  colnames(result) <- c("gene_symbol", "count")
  result$count <- as.integer(result$count)
  
  result
}

collapse_annotations <- function(data, key_col, value_col) {
  
  if (nrow(data) == 0) {
    return(data.frame(
      gene_symbol = character(),
      annotation = character()
    ))
  }
  
  valid <- !is.na(data[[key_col]]) &
    !is.na(data[[value_col]]) &
    nzchar(as.character(data[[value_col]]))
  
  data <- data[valid, , drop = FALSE]
  
  if (nrow(data) == 0) {
    return(data.frame(
      gene_symbol = character(),
      annotation = character()
    ))
  }
  
  split_values <- split(
    as.character(data[[value_col]]),
    as.character(data[[key_col]])
  )
  
  result <- data.frame(
    gene_symbol = names(split_values),
    annotation = vapply(
      split_values,
      function(x) paste(sort(unique(x)), collapse = "; "),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  
  result
}

safe_write_csv <- function(data, path) {
  
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  
  list_columns <- vapply(data, is.list, logical(1))
  
  if (any(list_columns)) {
    stop(
      "Cannot export list columns: ",
      paste(names(data)[list_columns], collapse = ", "),
      " in ",
      basename(path)
    )
  }
  
  write.csv(
    data,
    path,
    row.names = FALSE,
    na = ""
  )
}

# ------------------------------------------------------------
# 4. LOAD INPUTS
# ------------------------------------------------------------

cat("Loading input files...\n")

gene_pathway <- read_csv(
  required_files["gene_pathway_edges"]
)

gene_module <- read_csv(
  required_files["gene_module_edges"]
)

node_topology <- read_csv(
  required_files["node_topology"]
)

# ------------------------------------------------------------
# 5. IDENTIFY COLUMNS
# ------------------------------------------------------------

gp_gene_col <- find_column(
  gene_pathway,
  c("gene_symbol", "gene", "symbol"),
  TRUE,
  "gene_pathway_edges"
)

gp_pathway_col <- find_column(
  gene_pathway,
  c("pathway_id", "pathway", "term_id"),
  TRUE,
  "gene_pathway_edges"
)

gm_gene_col <- find_column(
  gene_module,
  c("gene_symbol", "gene", "symbol"),
  TRUE,
  "gene_module_edges"
)

gm_module_col <- find_column(
  gene_module,
  c("module", "module_name"),
  TRUE,
  "gene_module_edges"
)

topology_type_col <- find_column(
  node_topology,
  c("node_type"),
  TRUE,
  "node_topology_metrics"
)

topology_label_col <- find_column(
  node_topology,
  c("label", "node_label"),
  TRUE,
  "node_topology_metrics"
)

topology_betweenness_col <- find_column(
  node_topology,
  c("betweenness"),
  TRUE,
  "node_topology_metrics"
)

topology_degree_col <- find_column(
  node_topology,
  c("degree"),
  TRUE,
  "node_topology_metrics"
)

topology_articulation_col <- find_column(
  node_topology,
  c("is_articulation_point", "articulation_point"),
  FALSE,
  "node_topology_metrics"
)

cat("Input schemas validated.\n\n")

# ------------------------------------------------------------
# 6. NORMALIZE NETWORK TABLES
# ------------------------------------------------------------

gene_pathway_clean <- data.frame(
  gene_symbol = clean_gene(gene_pathway[[gp_gene_col]]),
  pathway_id = trimws(as.character(
    gene_pathway[[gp_pathway_col]]
  )),
  stringsAsFactors = FALSE
)

gene_module_clean <- data.frame(
  gene_symbol = clean_gene(gene_module[[gm_gene_col]]),
  module = trimws(as.character(
    gene_module[[gm_module_col]]
  )),
  stringsAsFactors = FALSE
)

gene_pathway_clean <- unique(
  gene_pathway_clean[
    !is.na(gene_pathway_clean$gene_symbol) &
      !is.na(gene_pathway_clean$pathway_id) &
      nzchar(gene_pathway_clean$pathway_id),
    ,
    drop = FALSE
  ]
)

gene_module_clean <- unique(
  gene_module_clean[
    !is.na(gene_module_clean$gene_symbol) &
      !is.na(gene_module_clean$module) &
      nzchar(gene_module_clean$module) &
      gene_module_clean$module != "NO_DEFINED_MODULE",
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 7. PATHWAY SUPPORT AND ANNOTATIONS
# ------------------------------------------------------------

pathway_support <- count_unique(
  gene_pathway_clean,
  "gene_symbol"
)

names(pathway_support)[2] <- "pathway_count"

pathway_annotations <- collapse_annotations(
  gene_pathway_clean,
  "gene_symbol",
  "pathway_id"
)

names(pathway_annotations)[2] <- "pathways"

# ------------------------------------------------------------
# 8. MODULE SUPPORT AND ANNOTATIONS
# ------------------------------------------------------------

module_support <- count_unique(
  gene_module_clean,
  "gene_symbol"
)

names(module_support)[2] <- "module_count"

module_annotations <- collapse_annotations(
  gene_module_clean,
  "gene_symbol",
  "module"
)

names(module_annotations)[2] <- "modules"

# ------------------------------------------------------------
# 9. EXTRACT GENE TOPOLOGY
# ------------------------------------------------------------

topology_type <- toupper(
  trimws(as.character(
    node_topology[[topology_type_col]]
  ))
)

topology_label <- clean_gene(
  node_topology[[topology_label_col]]
)

gene_rows <- topology_type == "GENE" &
  !is.na(topology_label)

gene_topology <- data.frame(
  gene_symbol = topology_label[gene_rows],
  degree = safe_numeric(
    node_topology[[topology_degree_col]][gene_rows]
  ),
  betweenness = safe_numeric(
    node_topology[[topology_betweenness_col]][gene_rows]
  ),
  stringsAsFactors = FALSE
)

if (!is.null(topology_articulation_col)) {
  
  articulation_values <- node_topology[
    gene_rows,
    topology_articulation_col
  ]
  
  gene_topology$is_articulation_point <- (
    tolower(as.character(articulation_values)) %in%
      c("true", "1", "yes")
  )
  
} else {
  
  gene_topology$is_articulation_point <- FALSE
  
}

gene_topology <- gene_topology[
  !is.na(gene_topology$gene_symbol),
  ,
  drop = FALSE
]

gene_topology <- gene_topology[
  !duplicated(gene_topology$gene_symbol),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 10. BRIDGE GENES
# ------------------------------------------------------------

bridge_genes <- character()

bridge_file <- file.path(
  topology_dir,
  "tables",
  "bridge_gene_topology.csv"
)

if (file.exists(bridge_file)) {
  
  bridge_data <- read_csv(bridge_file)
  
  bridge_gene_col <- find_column(
    bridge_data,
    c("gene_symbol", "gene", "label"),
    FALSE,
    "bridge_gene_topology"
  )
  
  if (!is.null(bridge_gene_col)) {
    bridge_genes <- unique(
      na.omit(clean_gene(bridge_data[[bridge_gene_col]]))
    )
  }
}

multi_module_genes <- module_support$gene_symbol[
  module_support$module_count > 1
]

bridge_genes <- unique(
  c(bridge_genes, multi_module_genes)
)

# ------------------------------------------------------------
# 11. BUILD CANDIDATE TABLE
# ------------------------------------------------------------

candidate_genes <- sort(unique(c(
  pathway_support$gene_symbol,
  module_support$gene_symbol,
  gene_topology$gene_symbol
)))

if (length(candidate_genes) == 0) {
  stop("No candidate genes could be extracted.")
}

candidates <- data.frame(
  gene_symbol = candidate_genes,
  stringsAsFactors = FALSE
)

for (source_table in list(
  pathway_support,
  module_support,
  gene_topology,
  pathway_annotations,
  module_annotations
)) {
  
  candidates <- merge(
    candidates,
    source_table,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )
}

# Explicitly convert annotation columns to character.
candidates$pathways <- as.character(candidates$pathways)
candidates$modules <- as.character(candidates$modules)

candidates$pathways[is.na(candidates$pathways)] <- ""
candidates$modules[is.na(candidates$modules)] <- ""

candidates$pathway_count[
  is.na(candidates$pathway_count)
] <- 0

candidates$module_count[
  is.na(candidates$module_count)
] <- 0

candidates$degree[
  !is.finite(candidates$degree)
] <- 0

candidates$betweenness[
  !is.finite(candidates$betweenness)
] <- 0

candidates$is_articulation_point[
  is.na(candidates$is_articulation_point)
] <- FALSE

candidates$is_bridge_gene <- (
  candidates$gene_symbol %in% bridge_genes
)

# ------------------------------------------------------------
# 12. NORMALIZE EVIDENCE
# ------------------------------------------------------------

candidates$betweenness_score <- scale_rank(
  candidates$betweenness
)

candidates$degree_score <- scale_rank(
  candidates$degree
)

candidates$pathway_score <- scale_rank(
  candidates$pathway_count
)

candidates$module_score <- scale_rank(
  candidates$module_count
)

candidates$bridge_score <- as.numeric(
  candidates$is_bridge_gene
)

candidates$articulation_score <- as.numeric(
  candidates$is_articulation_point
)

# ------------------------------------------------------------
# 13. COMPOSITE PRIORITY SCORE
# ------------------------------------------------------------

candidates$priority_score <- (
  0.30 * candidates$betweenness_score +
    0.20 * candidates$degree_score +
    0.15 * candidates$pathway_score +
    0.15 * candidates$module_score +
    0.10 * candidates$bridge_score +
    0.10 * candidates$articulation_score
) * 100

candidates$priority_class <- ifelse(
  candidates$priority_score >= 70,
  "HIGH",
  ifelse(
    candidates$priority_score >= 40,
    "MEDIUM",
    "EXPLORATORY"
  )
)

candidates <- candidates[
  order(
    -candidates$priority_score,
    -candidates$betweenness,
    -candidates$degree,
    candidates$gene_symbol
  ),
  ,
  drop = FALSE
]

rownames(candidates) <- NULL

candidates$priority_rank <- seq_len(nrow(candidates))

# ------------------------------------------------------------
# 14. EVIDENCE FLAGS
# ------------------------------------------------------------

candidates$network_supported <- (
  candidates$degree > 0 |
    candidates$betweenness > 0
)

candidates$multi_pathway_gene <- (
  candidates$pathway_count > 1
)

candidates$multi_module_gene <- (
  candidates$module_count > 1
)

candidates$high_betweenness <- (
  candidates$betweenness_score >= 0.90
)

candidates$high_degree <- (
  candidates$degree_score >= 0.90
)

# ------------------------------------------------------------
# 15. EXPORT COMPLETE RANKING
# ------------------------------------------------------------

safe_write_csv(
  candidates,
  file.path(table_dir, "candidate_gene_priority_ranking.csv")
)

# ------------------------------------------------------------
# 16. TOP 50 AND HIGH PRIORITY
# ------------------------------------------------------------

top_50 <- head(candidates, min(50, nrow(candidates)))

safe_write_csv(
  top_50,
  file.path(table_dir, "top_50_candidate_genes.csv")
)

high_priority <- candidates[
  candidates$priority_class == "HIGH",
  ,
  drop = FALSE
]

safe_write_csv(
  high_priority,
  file.path(table_dir, "high_priority_candidate_genes.csv")
)

# ------------------------------------------------------------
# 17. BRIDGE AND ARTICULATION RANKING
# ------------------------------------------------------------

bridge_ranking <- candidates[
  candidates$is_bridge_gene |
    candidates$is_articulation_point,
  ,
  drop = FALSE
]

bridge_ranking <- bridge_ranking[
  order(
    -bridge_ranking$priority_score,
    -bridge_ranking$betweenness
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  bridge_ranking,
  file.path(
    table_dir,
    "bridge_and_articulation_gene_ranking.csv"
  )
)

# ------------------------------------------------------------
# 18. MODULE PRIORITY SUMMARY
# ------------------------------------------------------------

module_gene_summary <- aggregate(
  gene_symbol ~ module,
  data = gene_module_clean,
  FUN = function(x) length(unique(x))
)

names(module_gene_summary)[2] <- "gene_count"

module_candidate_data <- merge(
  gene_module_clean,
  candidates[
    ,
    c("gene_symbol", "priority_score"),
    drop = FALSE
  ],
  by = "gene_symbol"
)

module_priority <- aggregate(
  priority_score ~ module,
  data = module_candidate_data,
  FUN = mean
)

module_priority <- merge(
  module_gene_summary,
  module_priority,
  by = "module",
  all = TRUE
)

module_priority <- module_priority[
  order(-module_priority$priority_score),
  ,
  drop = FALSE
]

safe_write_csv(
  module_priority,
  file.path(table_dir, "module_candidate_priority_summary.csv")
)

# ------------------------------------------------------------
# 19. PRIORITY CLASS SUMMARY
# ------------------------------------------------------------

priority_summary <- as.data.frame(
  table(
    factor(
      candidates$priority_class,
      levels = c("HIGH", "MEDIUM", "EXPLORATORY")
    )
  ),
  stringsAsFactors = FALSE
)

names(priority_summary) <- c(
  "priority_class",
  "gene_count"
)

priority_summary$percentage <- (
  priority_summary$gene_count /
    sum(priority_summary$gene_count)
) * 100

safe_write_csv(
  priority_summary,
  file.path(table_dir, "priority_class_summary.csv")
)

# ------------------------------------------------------------
# 20. EVIDENCE MATRIX
# ------------------------------------------------------------

evidence_columns <- c(
  "priority_rank",
  "gene_symbol",
  "priority_score",
  "priority_class",
  "degree",
  "betweenness",
  "pathway_count",
  "module_count",
  "is_bridge_gene",
  "is_articulation_point",
  "pathways",
  "modules"
)

evidence_table <- candidates[
  ,
  evidence_columns,
  drop = FALSE
]

safe_write_csv(
  evidence_table,
  file.path(table_dir, "candidate_gene_evidence_matrix.csv")
)

# ------------------------------------------------------------
# 21. REPORT
# ------------------------------------------------------------

top10 <- head(candidates, 10)

top10_text <- paste(
  sprintf(
    "%d. %s | Score: %.2f | %s | Degree: %.0f | Betweenness: %.5f | Pathways: %.0f | Modules: %.0f",
    top10$priority_rank,
    top10$gene_symbol,
    top10$priority_score,
    top10$priority_class,
    top10$degree,
    top10$betweenness,
    top10$pathway_count,
    top10$module_count
  ),
  collapse = "\n"
)

bridge_text <- if (nrow(bridge_ranking) > 0) {
  paste(
    sprintf(
      "%s | Score: %.2f | Betweenness: %.5f | Modules: %.0f",
      bridge_ranking$gene_symbol,
      bridge_ranking$priority_score,
      bridge_ranking$betweenness,
      bridge_ranking$module_count
    ),
    collapse = "\n"
  )
} else {
  "No bridge or articulation genes identified."
}

report <- c(
  "TRANSCRIPTOME-X",
  "NETWORK-BASED CANDIDATE GENE PRIORITIZATION",
  "============================================================",
  "",
  paste("Total candidate genes:", nrow(candidates)),
  paste("High priority:", sum(candidates$priority_class == "HIGH")),
  paste("Medium priority:", sum(candidates$priority_class == "MEDIUM")),
  paste("Exploratory:", sum(candidates$priority_class == "EXPLORATORY")),
  paste("Bridge genes:", sum(candidates$is_bridge_gene)),
  paste("Articulation genes:", sum(candidates$is_articulation_point)),
  "",
  "TOP 10 CANDIDATE GENES",
  "----------------------",
  top10_text,
  "",
  "BRIDGE AND ARTICULATION GENES",
  "-----------------------------",
  bridge_text,
  "",
  "SCORING FRAMEWORK",
  "-----------------",
  "Betweenness centrality: 30%",
  "Degree: 20%",
  "Pathway support: 15%",
  "Module support: 15%",
  "Bridge gene status: 10%",
  "Articulation point status: 10%",
  "",
  "INTERPRETATION",
  "--------------",
  "The priority score is a descriptive composite ranking.",
  "It combines network topology and annotation coverage.",
  "It is not a statistical p-value or a measure of causal relevance.",
  "Genes with more annotations may receive higher support scores.",
  "The ranking should be interpreted alongside differential",
  "expression, pathway enrichment and experimental evidence.",
  "",
  "LIMITATIONS",
  "-----------",
  "The network is based on gene-pathway-module annotations.",
  "Edges are not necessarily physical molecular interactions.",
  "Pathway and module counts depend on selected annotations.",
  "Centrality depends on the constructed network topology.",
  "This analysis does not demonstrate pathway activation."
)

writeLines(
  report,
  file.path(report_dir, "candidate_gene_prioritization_report.txt")
)

# ------------------------------------------------------------
# 22. AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "input_gene_pathway_records",
    "unique_gene_pathway_records",
    "input_gene_module_records",
    "unique_gene_module_records",
    "topology_gene_records",
    "candidate_genes",
    "genes_with_pathway_support",
    "genes_with_module_support",
    "genes_with_topology_support",
    "bridge_genes",
    "articulation_genes",
    "high_priority_genes",
    "medium_priority_genes",
    "exploratory_genes"
  ),
  value = c(
    nrow(gene_pathway),
    nrow(gene_pathway_clean),
    nrow(gene_module),
    nrow(gene_module_clean),
    nrow(gene_topology),
    nrow(candidates),
    sum(candidates$pathway_count > 0),
    sum(candidates$module_count > 0),
    sum(candidates$degree > 0 | candidates$betweenness > 0),
    sum(candidates$is_bridge_gene),
    sum(candidates$is_articulation_point),
    sum(candidates$priority_class == "HIGH"),
    sum(candidates$priority_class == "MEDIUM"),
    sum(candidates$priority_class == "EXPLORATORY")
  )
)

safe_write_csv(
  audit,
  file.path(audit_dir, "script40_audit.csv")
)

# ------------------------------------------------------------
# 23. FINAL VALIDATION
# ------------------------------------------------------------

if (anyDuplicated(candidates$gene_symbol)) {
  stop("Validation failed: duplicate candidate genes.")
}

if (any(!is.finite(candidates$priority_score))) {
  stop("Validation failed: non-finite priority scores.")
}

if (any(
  candidates$priority_score < 0 |
  candidates$priority_score > 100
)) {
  stop("Validation failed: priority scores outside 0-100.")
}

if (!identical(
  candidates$priority_rank,
  seq_len(nrow(candidates))
)) {
  stop("Validation failed: incorrect ranking sequence.")
}

all_output_files <- list.files(
  output_dir,
  recursive = TRUE,
  full.names = TRUE
)

cat("\n============================================================\n")
cat("SCRIPT 40 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Candidate genes:", nrow(candidates), "\n")
cat("High priority:", sum(candidates$priority_class == "HIGH"), "\n")
cat("Medium priority:", sum(candidates$priority_class == "MEDIUM"), "\n")
cat("Exploratory:", sum(candidates$priority_class == "EXPLORATORY"), "\n")
cat("Bridge genes:", sum(candidates$is_bridge_gene), "\n")
cat("Articulation genes:", sum(candidates$is_articulation_point), "\n")

cat("\nTOP 10 CANDIDATES:\n")

print(
  top10[
    ,
    c(
      "priority_rank",
      "gene_symbol",
      "priority_score",
      "priority_class"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nOutput files:", length(all_output_files), "\n")
cat("Output directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
