
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 37 - GENE-PATHWAY-MODULE NETWORK ANALYSIS
# CORRECTED VERSION
#
# Objectives:
# 1. Build gene-pathway-module network tables.
# 2. Correctly expand multi-module gene assignments.
# 3. Identify multi-module bridge genes.
# 4. Rank gene hub candidates.
# 5. Summarize network structure.
#
# INPUT:
# results/pathway_concordance_specificity/tables/
#   pathway_gene_evidence.csv
#
# results/shared_pathway_priority/tables/
#   pathway_priority_ranking.csv
#
# OUTPUT:
# results/gene_pathway_module_network/
#
# Important:
# NO_DEFINED_MODULE is not a biological module.
# Module assignments separated by "|" are expanded.
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

PROJECT_ROOT <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = FALSE
)

EVIDENCE_FILE <- file.path(
  PROJECT_ROOT,
  "results",
  "pathway_concordance_specificity",
  "tables",
  "pathway_gene_evidence.csv"
)

PRIORITY_FILE <- file.path(
  PROJECT_ROOT,
  "results",
  "shared_pathway_priority",
  "tables",
  "pathway_priority_ranking.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_pathway_module_network"
)

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
AUDIT_DIR <- file.path(OUTPUT_DIR, "audit")
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")

for (directory in c(
  TABLE_DIR,
  AUDIT_DIR,
  REPORT_DIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!file.exists(EVIDENCE_FILE)) {
  stop("Missing input: ", EVIDENCE_FILE)
}

if (!file.exists(PRIORITY_FILE)) {
  stop("Missing input: ", PRIORITY_FILE)
}


# ============================================================
# 2. UTILITY FUNCTIONS
# ============================================================

safe_character <- function(x) {
  
  result <- trimws(as.character(x))
  result[is.na(result)] <- ""
  
  result
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(as.character(x))
  )
}


safe_logical <- function(x) {
  
  if (is.logical(x)) {
    return(ifelse(is.na(x), FALSE, x))
  }
  
  values <- toupper(
    trimws(as.character(x))
  )
  
  values[is.na(values)] <- ""
  
  values %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )
}


write_csv <- function(data, filename) {
  
  write.csv(
    data,
    filename,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}


split_modules <- function(x) {
  
  x <- safe_character(x)
  
  if (length(x) == 0 || !nzchar(x[1])) {
    return(character(0))
  }
  
  modules <- unlist(
    strsplit(
      x,
      "\\s*\\|\\s*",
      perl = TRUE
    ),
    use.names = FALSE
  )
  
  modules <- trimws(modules)
  
  modules <- unique(
    modules[
      nzchar(modules) &
        modules != "NO_DEFINED_MODULE"
    ]
  )
  
  sort(modules)
}


collapse_values <- function(x) {
  
  values <- sort(
    unique(
      safe_character(x)
    )
  )
  
  values <- values[nzchar(values)]
  
  if (length(values) == 0) {
    return("")
  }
  
  paste(values, collapse = ";")
}


minmax_score <- function(x) {
  
  x <- safe_numeric(x)
  
  result <- rep(0, length(x))
  
  valid <- is.finite(x)
  
  if (!any(valid)) {
    return(result)
  }
  
  xmin <- min(x[valid])
  xmax <- max(x[valid])
  
  if (xmax == xmin) {
    result[valid] <- 100
  } else {
    result[valid] <- (
      (x[valid] - xmin) /
        (xmax - xmin)
    ) * 100
  }
  
  result[!is.finite(result)] <- 0
  
  result
}


# ============================================================
# 3. LOAD INPUTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 37 - GENE-PATHWAY-MODULE NETWORK ANALYSIS\n")
cat("CORRECTED VERSION\n")
cat("============================================================\n\n")

cat("Loading pathway-gene evidence...\n")

evidence <- read.csv(
  EVIDENCE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

cat("Loading pathway priority ranking...\n")

priority <- read.csv(
  PRIORITY_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_evidence <- c(
  "gene_symbol",
  "pathway_id",
  "term",
  "ontology",
  "group",
  "modules",
  "significant_siNRP1_3",
  "significant_siNRP1_5",
  "significant_both",
  "direction_siNRP1_3",
  "direction_siNRP1_5"
)

missing_evidence <- setdiff(
  required_evidence,
  colnames(evidence)
)

if (length(missing_evidence) > 0) {
  stop(
    "Missing evidence columns: ",
    paste(missing_evidence, collapse = ", ")
  )
}

required_priority <- c(
  "pathway_id",
  "pathway_status",
  "priority_score",
  "priority_rank"
)

missing_priority <- setdiff(
  required_priority,
  colnames(priority)
)

if (length(missing_priority) > 0) {
  stop(
    "Missing priority columns: ",
    paste(missing_priority, collapse = ", ")
  )
}


# ============================================================
# 4. NORMALIZE DATA
# ============================================================

evidence$gene_symbol <- safe_character(
  evidence$gene_symbol
)

evidence$pathway_id <- safe_character(
  evidence$pathway_id
)

evidence$term <- safe_character(
  evidence$term
)

evidence$ontology <- safe_character(
  evidence$ontology
)

evidence$group <- safe_character(
  evidence$group
)

evidence$modules <- safe_character(
  evidence$modules
)

evidence$significant_siNRP1_3 <- safe_logical(
  evidence$significant_siNRP1_3
)

evidence$significant_siNRP1_5 <- safe_logical(
  evidence$significant_siNRP1_5
)

evidence$significant_both <- safe_logical(
  evidence$significant_both
)

evidence$direction_siNRP1_3 <- toupper(
  safe_character(evidence$direction_siNRP1_3)
)

evidence$direction_siNRP1_5 <- toupper(
  safe_character(evidence$direction_siNRP1_5)
)

priority$pathway_id <- safe_character(
  priority$pathway_id
)

priority$pathway_status <- safe_character(
  priority$pathway_status
)

priority$priority_score <- safe_numeric(
  priority$priority_score
)

priority$priority_rank <- safe_numeric(
  priority$priority_rank
)

evidence <- evidence[
  nzchar(evidence$gene_symbol) &
    nzchar(evidence$pathway_id),
  ,
  drop = FALSE
]

priority <- priority[
  nzchar(priority$pathway_id),
  ,
  drop = FALSE
]

if (nrow(evidence) == 0) {
  stop("No valid evidence records.")
}

if (anyDuplicated(priority$pathway_id)) {
  stop("Priority table contains duplicated pathway_id values.")
}


# ============================================================
# 5. REMOVE DUPLICATE GENE-PATHWAY RECORDS
# ============================================================

duplicate_count <- sum(
  duplicated(
    evidence[, c(
      "gene_symbol",
      "pathway_id"
    )]
  )
)

evidence <- evidence[
  !duplicated(
    evidence[, c(
      "gene_symbol",
      "pathway_id"
    )]
  ),
  ,
  drop = FALSE
]

cat("Evidence records:", nrow(evidence), "\n")
cat(
  "Unique genes:",
  length(unique(evidence$gene_symbol)),
  "\n"
)
cat(
  "Unique pathways:",
  length(unique(evidence$pathway_id)),
  "\n"
)
cat(
  "Duplicate records removed:",
  duplicate_count,
  "\n\n"
)


# ============================================================
# 6. EXPAND GENE-MODULE ASSIGNMENTS
# ============================================================

module_records <- list()

for (i in seq_len(nrow(evidence))) {
  
  modules <- split_modules(
    evidence$modules[i]
  )
  
  if (length(modules) == 0) {
    next
  }
  
  for (module in modules) {
    
    module_records[[length(module_records) + 1]] <-
      data.frame(
        gene_symbol = evidence$gene_symbol[i],
        module = module,
        pathway_id = evidence$pathway_id[i],
        significant_siNRP1_3 =
          evidence$significant_siNRP1_3[i],
        significant_siNRP1_5 =
          evidence$significant_siNRP1_5[i],
        significant_both =
          evidence$significant_both[i],
        direction_siNRP1_3 =
          evidence$direction_siNRP1_3[i],
        direction_siNRP1_5 =
          evidence$direction_siNRP1_5[i],
        stringsAsFactors = FALSE
      )
  }
}

if (length(module_records) > 0) {
  
  gene_module_edges <- do.call(
    rbind,
    module_records
  )
  
  gene_module_edges <- gene_module_edges[
    !duplicated(
      gene_module_edges[, c(
        "gene_symbol",
        "module"
      )]
    ),
    ,
    drop = FALSE
  ]
  
  rownames(gene_module_edges) <- NULL
  
} else {
  
  gene_module_edges <- data.frame(
    gene_symbol = character(0),
    module = character(0),
    pathway_id = character(0),
    significant_siNRP1_3 = logical(0),
    significant_siNRP1_5 = logical(0),
    significant_both = logical(0),
    direction_siNRP1_3 = character(0),
    direction_siNRP1_5 = character(0),
    stringsAsFactors = FALSE
  )
}

cat(
  "Expanded gene-module associations:",
  nrow(gene_module_edges),
  "\n"
)

cat(
  "Unique biological modules:",
  length(unique(gene_module_edges$module)),
  "\n\n"
)


# ============================================================
# 7. GENE-PATHWAY EDGES
# ============================================================

gene_pathway_edges <- merge(
  evidence[, c(
    "gene_symbol",
    "pathway_id",
    "term",
    "ontology",
    "group",
    "significant_siNRP1_3",
    "significant_siNRP1_5",
    "significant_both",
    "direction_siNRP1_3",
    "direction_siNRP1_5"
  )],
  priority[, c(
    "pathway_id",
    "pathway_status",
    "priority_score",
    "priority_rank"
  )],
  by = "pathway_id",
  all.x = TRUE,
  sort = FALSE
)

gene_pathway_edges <- gene_pathway_edges[
  order(
    gene_pathway_edges$gene_symbol,
    gene_pathway_edges$pathway_id
  ),
  ,
  drop = FALSE
]

rownames(gene_pathway_edges) <- NULL


# ============================================================
# 8. PATHWAY-MODULE EDGES
# ============================================================

pathway_module_edges <- unique(
  gene_module_edges[, c(
    "pathway_id",
    "module"
  )]
)

pathway_module_edges <- pathway_module_edges[
  order(
    pathway_module_edges$pathway_id,
    pathway_module_edges$module
  ),
  ,
  drop = FALSE
]

rownames(pathway_module_edges) <- NULL


# ============================================================
# 9. GENE-LEVEL NETWORK METRICS
# ============================================================

genes <- sort(
  unique(evidence$gene_symbol)
)

gene_records <- lapply(
  genes,
  function(gene) {
    
    d <- evidence[
      evidence$gene_symbol == gene,
      ,
      drop = FALSE
    ]
    
    module_d <- gene_module_edges[
      gene_module_edges$gene_symbol == gene,
      ,
      drop = FALSE
    ]
    
    gene_modules <- sort(
      unique(module_d$module)
    )
    
    pathway_status <- merge(
      d[, "pathway_id", drop = FALSE],
      priority[, c(
        "pathway_id",
        "pathway_status",
        "priority_score"
      )],
      by = "pathway_id",
      all.x = TRUE
    )
    
    sig3 <- any(d$significant_siNRP1_3)
    sig5 <- any(d$significant_siNRP1_5)
    sigboth <- any(d$significant_both)
    
    directions <- d[
      d$significant_both &
        d$direction_siNRP1_3 %in% c("UP", "DOWN") &
        d$direction_siNRP1_5 %in% c("UP", "DOWN"),
      ,
      drop = FALSE
    ]
    
    concordance <- if (nrow(directions) == 0) {
      "NOT_EVALUABLE"
    } else if (all(
      directions$direction_siNRP1_3 ==
      directions$direction_siNRP1_5
    )) {
      "CONCORDANT"
    } else if (all(
      directions$direction_siNRP1_3 !=
      directions$direction_siNRP1_5
    )) {
      "DISCORDANT"
    } else {
      "MIXED"
    }
    
    data.frame(
      gene_symbol = gene,
      pathway_count = length(unique(d$pathway_id)),
      module_count = length(gene_modules),
      shared_pathway_count = sum(
        pathway_status$pathway_status == "SHARED",
        na.rm = TRUE
      ),
      siNRP1_3_specific_pathway_count = sum(
        pathway_status$pathway_status == "ONLY_siNRP1_3",
        na.rm = TRUE
      ),
      siNRP1_5_specific_pathway_count = sum(
        pathway_status$pathway_status == "ONLY_siNRP1_5",
        na.rm = TRUE
      ),
      significant_siNRP1_3 = sig3,
      significant_siNRP1_5 = sig5,
      significant_both = sigboth,
      direction_siNRP1_3 = collapse_values(
        d$direction_siNRP1_3[
          d$significant_siNRP1_3
        ]
      ),
      direction_siNRP1_5 = collapse_values(
        d$direction_siNRP1_5[
          d$significant_siNRP1_5
        ]
      ),
      concordance_status = concordance,
      pathway_list = collapse_values(
        paste(
          d$ontology,
          d$group,
          d$term,
          sep = "|||"
        )
      ),
      module_list = paste(
        gene_modules,
        collapse = " | "
      ),
      mean_pathway_priority = if (
        all(!is.finite(pathway_status$priority_score))
      ) {
        NA_real_
      } else {
        mean(
          pathway_status$priority_score,
          na.rm = TRUE
        )
      },
      maximum_pathway_priority = if (
        all(!is.finite(pathway_status$priority_score))
      ) {
        NA_real_
      } else {
        max(
          pathway_status$priority_score,
          na.rm = TRUE
        )
      },
      stringsAsFactors = FALSE
    )
  }
)

gene_metrics <- do.call(
  rbind,
  gene_records
)

rownames(gene_metrics) <- NULL


# ============================================================
# 10. HUB SCORES
# ============================================================

gene_metrics$score_pathway_recurrence <- minmax_score(
  log1p(gene_metrics$pathway_count)
)

gene_metrics$score_module_diversity <- minmax_score(
  gene_metrics$module_count
)

gene_metrics$score_significance <- (
  as.numeric(gene_metrics$significant_siNRP1_3) +
    as.numeric(gene_metrics$significant_siNRP1_5) +
    as.numeric(gene_metrics$significant_both)
) / 3 * 100

gene_metrics$score_shared_pathways <- minmax_score(
  gene_metrics$shared_pathway_count
)

gene_metrics$score_pathway_priority <- minmax_score(
  gene_metrics$mean_pathway_priority
)

gene_metrics$hub_score <-
  0.30 * gene_metrics$score_pathway_recurrence +
  0.20 * gene_metrics$score_module_diversity +
  0.20 * gene_metrics$score_significance +
  0.15 * gene_metrics$score_shared_pathways +
  0.15 * gene_metrics$score_pathway_priority

gene_metrics$hub_score <- round(
  gene_metrics$hub_score,
  2
)

gene_metrics <- gene_metrics[
  order(
    -gene_metrics$hub_score,
    -gene_metrics$pathway_count,
    gene_metrics$gene_symbol
  ),
  ,
  drop = FALSE
]

gene_metrics$hub_rank <- seq_len(
  nrow(gene_metrics)
)

gene_metrics$hub_class <- ifelse(
  gene_metrics$hub_rank <= ceiling(
    nrow(gene_metrics) * 0.10
  ),
  "TOP_10_PERCENT",
  ifelse(
    gene_metrics$hub_rank <= ceiling(
      nrow(gene_metrics) * 0.25
    ),
    "TOP_25_PERCENT",
    "OTHER"
  )
)

gene_metrics$bridge_status <- ifelse(
  gene_metrics$module_count > 1,
  "MULTI_MODULE_BRIDGE",
  ifelse(
    gene_metrics$module_count == 1,
    "SINGLE_MODULE",
    "NO_DEFINED_MODULE"
  )
)

gene_metrics$cross_contrast_status <-
  gene_metrics$concordance_status


# ============================================================
# 11. BRIDGE GENE TABLE
# ============================================================

bridge_genes <- gene_metrics[
  gene_metrics$module_count > 1,
  ,
  drop = FALSE
]

bridge_genes <- bridge_genes[
  order(
    -bridge_genes$module_count,
    -bridge_genes$hub_score
  ),
  ,
  drop = FALSE
]

rownames(bridge_genes) <- NULL


# ============================================================
# 12. SHARED PATHWAY GENE CANDIDATES
# ============================================================

shared_candidates <- gene_metrics[
  gene_metrics$shared_pathway_count > 0,
  ,
  drop = FALSE
]

shared_candidates <- shared_candidates[
  order(
    -shared_candidates$shared_pathway_count,
    -shared_candidates$hub_score
  ),
  ,
  drop = FALSE
]

rownames(shared_candidates) <- NULL


# ============================================================
# 13. MODULE SUMMARY
# ============================================================

modules <- sort(
  unique(gene_module_edges$module)
)

module_records <- lapply(
  modules,
  function(module_name) {
    
    d <- gene_module_edges[
      gene_module_edges$module == module_name,
      ,
      drop = FALSE
    ]
    
    module_genes <- unique(
      d$gene_symbol
    )
    
    metrics <- gene_metrics[
      gene_metrics$gene_symbol %in% module_genes,
      ,
      drop = FALSE
    ]
    
    pathway_count <- length(
      unique(d$pathway_id)
    )
    
    sig_both <- unique(
      d$gene_symbol[d$significant_both]
    )
    
    concordant <- metrics$gene_symbol[
      metrics$concordance_status == "CONCORDANT"
    ]
    
    discordant <- metrics$gene_symbol[
      metrics$concordance_status == "DISCORDANT"
    ]
    
    data.frame(
      module = module_name,
      unique_genes = length(module_genes),
      pathways_represented = pathway_count,
      significant_both = length(sig_both),
      concordant_genes = sum(
        module_genes %in% concordant
      ),
      discordant_genes = sum(
        module_genes %in% discordant
      ),
      multi_module_genes = sum(
        metrics$module_count > 1
      ),
      mean_hub_score = round(
        mean(metrics$hub_score),
        2
      ),
      stringsAsFactors = FALSE
    )
  }
)

network_module_summary <- if (
  length(module_records) > 0
) {
  do.call(rbind, module_records)
} else {
  data.frame(
    module = character(0),
    unique_genes = integer(0),
    pathways_represented = integer(0),
    significant_both = integer(0),
    concordant_genes = integer(0),
    discordant_genes = integer(0),
    multi_module_genes = integer(0),
    mean_hub_score = numeric(0)
  )
}

network_module_summary <- network_module_summary[
  order(
    -network_module_summary$unique_genes
  ),
  ,
  drop = FALSE
]


# ============================================================
# 14. NETWORK SUMMARY
# ============================================================

network_summary <- data.frame(
  metric = c(
    "gene_nodes",
    "pathway_nodes",
    "module_nodes",
    "gene_pathway_edges",
    "gene_module_edges",
    "pathway_module_edges",
    "multi_module_bridge_genes",
    "hub_candidates_top_10_percent",
    "genes_without_defined_module"
  ),
  value = c(
    length(unique(evidence$gene_symbol)),
    length(unique(evidence$pathway_id)),
    length(modules),
    nrow(gene_pathway_edges),
    nrow(gene_module_edges),
    nrow(pathway_module_edges),
    nrow(bridge_genes),
    ceiling(nrow(gene_metrics) * 0.10),
    sum(gene_metrics$module_count == 0)
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 15. VALIDATION OF KNOWN BRIDGE GENES
# ============================================================

known_bridge_genes <- c(
  "CCN2",
  "CCN3",
  "THBS1"
)

known_bridge_validation <- gene_metrics[
  gene_metrics$gene_symbol %in% known_bridge_genes,
  c(
    "gene_symbol",
    "module_count",
    "module_list",
    "bridge_status",
    "hub_score",
    "hub_rank"
  ),
  drop = FALSE
]

expected_counts <- c(
  CCN2 = 3,
  CCN3 = 2,
  THBS1 = 2
)

if (nrow(known_bridge_validation) != 3) {
  stop(
    "Validation failure: one or more known bridge genes are missing."
  )
}

for (i in seq_len(nrow(known_bridge_validation))) {
  
  gene <- known_bridge_validation$gene_symbol[i]
  
  expected <- expected_counts[[gene]]
  
  observed <- known_bridge_validation$module_count[i]
  
  if (observed != expected) {
    stop(
      "Bridge validation failed for ",
      gene,
      ": expected ",
      expected,
      " modules, observed ",
      observed
    )
  }
}


# ============================================================
# 16. AUDIT
# ============================================================

audit <- data.frame(
  metric = c(
    "input_evidence_records",
    "duplicate_records_removed",
    "unique_genes",
    "unique_pathways",
    "unique_biological_modules",
    "expanded_gene_module_edges",
    "gene_pathway_edges",
    "pathway_module_edges",
    "multi_module_bridge_genes",
    "genes_without_defined_module",
    "hub_candidates_top_10_percent",
    "CCN2_module_count",
    "CCN3_module_count",
    "THBS1_module_count",
    "known_bridge_validation"
  ),
  value = c(
    nrow(evidence),
    duplicate_count,
    length(unique(evidence$gene_symbol)),
    length(unique(evidence$pathway_id)),
    length(modules),
    nrow(gene_module_edges),
    nrow(gene_pathway_edges),
    nrow(pathway_module_edges),
    nrow(bridge_genes),
    sum(gene_metrics$module_count == 0),
    ceiling(nrow(gene_metrics) * 0.10),
    expected_counts[["CCN2"]],
    expected_counts[["CCN3"]],
    expected_counts[["THBS1"]],
    "PASSED"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 17. EXPORT TABLES
# ============================================================

write_csv(
  gene_metrics,
  file.path(
    TABLE_DIR,
    "gene_hub_ranking.csv"
  )
)

write_csv(
  bridge_genes,
  file.path(
    TABLE_DIR,
    "multi_module_bridge_genes.csv"
  )
)

write_csv(
  shared_candidates,
  file.path(
    TABLE_DIR,
    "shared_pathway_gene_candidates.csv"
  )
)

write_csv(
  gene_pathway_edges,
  file.path(
    TABLE_DIR,
    "gene_pathway_edges.csv"
  )
)

write_csv(
  gene_module_edges,
  file.path(
    TABLE_DIR,
    "gene_module_edges.csv"
  )
)

write_csv(
  pathway_module_edges,
  file.path(
    TABLE_DIR,
    "pathway_module_edges.csv"
  )
)

write_csv(
  network_module_summary,
  file.path(
    TABLE_DIR,
    "network_module_summary.csv"
  )
)

write_csv(
  network_summary,
  file.path(
    TABLE_DIR,
    "network_summary.csv"
  )
)

write_csv(
  known_bridge_validation,
  file.path(
    TABLE_DIR,
    "known_bridge_gene_validation.csv"
  )
)

write_csv(
  audit,
  file.path(
    AUDIT_DIR,
    "script37_audit.csv"
  )
)


# ============================================================
# 18. REPORT
# ============================================================

report_file <- file.path(
  REPORT_DIR,
  "gene_pathway_module_network_report.txt"
)

sink(report_file)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 37 - GENE-PATHWAY-MODULE NETWORK ANALYSIS\n")
cat("CORRECTED VERSION\n")
cat("============================================================\n\n")

cat("NETWORK SIZE\n")
print(network_summary, row.names = FALSE)

cat("\nMODULE SUMMARY\n")
print(network_module_summary, row.names = FALSE)

cat("\nKNOWN BRIDGE GENE VALIDATION\n")
print(known_bridge_validation, row.names = FALSE)

cat("\nTOP 20 HUB CANDIDATES\n")

print(
  head(
    gene_metrics[, c(
      "gene_symbol",
      "pathway_count",
      "module_count",
      "shared_pathway_count",
      "hub_score",
      "hub_rank",
      "hub_class",
      "bridge_status"
    )],
    20
  ),
  row.names = FALSE
)

cat("\nMULTI-MODULE BRIDGE GENES\n")

if (nrow(bridge_genes) > 0) {
  
  print(
    bridge_genes[, c(
      "gene_symbol",
      "module_count",
      "module_list",
      "pathway_count",
      "hub_score",
      "hub_rank"
    )],
    row.names = FALSE
  )
  
} else {
  
  cat("No multi-module bridge genes identified.\n")
}

cat("\nMETHODOLOGICAL NOTES\n")
cat(
  "- Multi-module annotations are expanded into separate associations.\n"
)
cat(
  "- NO_DEFINED_MODULE is excluded from biological module counts.\n"
)
cat(
  "- Hub scores are descriptive ranking metrics.\n"
)
cat(
  "- Gene-pathway and gene-module edges represent annotations, not causal interactions.\n"
)
cat(
  "- A gene may be associated with multiple pathways and modules.\n"
)
cat(
  "- No new differential expression or enrichment analysis was performed.\n"
)

cat("\n============================================================\n")
cat("END OF SCRIPT 37\n")
cat("============================================================\n")

sink()


# ============================================================
# 19. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 37 COMPLETED - CORRECTED VERSION\n")
cat("============================================================\n\n")

cat("Gene nodes:", length(unique(evidence$gene_symbol)), "\n")
cat("Pathway nodes:", length(unique(evidence$pathway_id)), "\n")
cat("Module nodes:", length(modules), "\n")
cat("Gene-pathway edges:", nrow(gene_pathway_edges), "\n")
cat("Gene-module edges:", nrow(gene_module_edges), "\n")
cat("Pathway-module edges:", nrow(pathway_module_edges), "\n")
cat("Multi-module bridge genes:", nrow(bridge_genes), "\n")
cat(
  "Top 10% hub candidates:",
  ceiling(nrow(gene_metrics) * 0.10),
  "\n"
)

cat("\nKNOWN BRIDGE GENES\n")
print(
  known_bridge_validation,
  row.names = FALSE
)

cat("\nOUTPUT FILES\n")
cat(
  "1. gene_hub_ranking.csv\n",
  "2. multi_module_bridge_genes.csv\n",
  "3. shared_pathway_gene_candidates.csv\n",
  "4. gene_pathway_edges.csv\n",
  "5. gene_module_edges.csv\n",
  "6. pathway_module_edges.csv\n",
  "7. network_module_summary.csv\n",
  "8. network_summary.csv\n",
  "9. known_bridge_gene_validation.csv\n",
  "10. script37_audit.csv\n",
  "11. gene_pathway_module_network_report.txt\n",
  sep = ""
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n============================================================\n")
