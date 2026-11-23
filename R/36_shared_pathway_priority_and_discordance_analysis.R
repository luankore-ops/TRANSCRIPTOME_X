
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 36
# SHARED PATHWAY PRIORITY AND DISCORDANCE ANALYSIS
#
# Objective:
# 1. Prioritize pathways identified in Script 35.
# 2. Characterize shared pathways.
# 3. Quantify gene-level concordance and discordance.
# 4. Identify pathways with stronger evidence in each knockdown.
# 5. Summarize pathway-associated functional modules.
#
# INPUT:
# results/pathway_concordance_specificity/
#
# OUTPUT:
# results/shared_pathway_priority/
#
# Methodological note:
# This is a descriptive prioritization analysis.
# It does not perform new enrichment tests and does not
# establish pathway activation, inhibition, or causality.
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

PROJECT_ROOT <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = FALSE
)

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "pathway_concordance_specificity"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "shared_pathway_priority"
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

INPUT_FILE <- file.path(
  INPUT_DIR,
  "tables",
  "pathway_gene_evidence.csv"
)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found: ",
    INPUT_FILE
  )
}


# ============================================================
# 2. UTILITIES
# ============================================================

write_csv <- function(data, filename) {
  
  write.csv(
    data,
    file = filename,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
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


safe_character <- function(x) {
  
  result <- trimws(as.character(x))
  
  result[is.na(result)] <- ""
  
  result
}


unique_nonempty <- function(x) {
  
  x <- safe_character(x)
  
  sort(
    unique(
      x[nzchar(x)]
    )
  )
}


collapse_unique <- function(x) {
  
  values <- unique_nonempty(x)
  
  if (length(values) == 0) {
    return("")
  }
  
  paste(values, collapse = ";")
}


parse_modules <- function(x) {
  
  x <- safe_character(x)
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  result <- unlist(
    strsplit(x, ";", fixed = TRUE),
    use.names = FALSE
  )
  
  unique_nonempty(result)
}


minmax_score <- function(x) {
  
  x <- safe_numeric(x)
  
  result <- rep(0, length(x))
  
  valid <- is.finite(x)
  
  if (!any(valid)) {
    return(result)
  }
  
  minimum <- min(x[valid])
  maximum <- max(x[valid])
  
  if (maximum == minimum) {
    result[valid] <- 100
  } else {
    result[valid] <- (
      (x[valid] - minimum) /
        (maximum - minimum)
    ) * 100
  }
  
  result[!is.finite(result)] <- 0
  
  result
}


safe_fraction <- function(numerator, denominator) {
  
  result <- rep(0, length(numerator))
  
  valid <- denominator > 0
  
  result[valid] <-
    numerator[valid] / denominator[valid]
  
  result
}


# ============================================================
# 3. LOAD AND VALIDATE INPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 36 - PATHWAY PRIORITY ANALYSIS\n")
cat("============================================================\n\n")

cat("Input:", INPUT_FILE, "\n")

evidence <- read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_columns <- c(
  "gene_symbol",
  "pathway_id",
  "term",
  "ontology",
  "group",
  "logFC_siNRP1_3",
  "FDR_siNRP1_3",
  "significant_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_5",
  "significant_siNRP1_5",
  "in_signature_siNRP1_3",
  "in_signature_siNRP1_5",
  "significant_both",
  "direction_siNRP1_3",
  "direction_siNRP1_5",
  "joint_class",
  "modules"
)

missing_columns <- setdiff(
  required_columns,
  colnames(evidence)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(evidence) == 0) {
  stop("The input evidence file contains no rows.")
}

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

evidence <- evidence[
  nzchar(evidence$gene_symbol) &
    nzchar(evidence$pathway_id),
  ,
  drop = FALSE
]

if (nrow(evidence) == 0) {
  stop("No valid pathway-gene records remain after filtering.")
}

evidence$logFC_siNRP1_3 <- safe_numeric(
  evidence$logFC_siNRP1_3
)

evidence$FDR_siNRP1_3 <- safe_numeric(
  evidence$FDR_siNRP1_3
)

evidence$logFC_siNRP1_5 <- safe_numeric(
  evidence$logFC_siNRP1_5
)

evidence$FDR_siNRP1_5 <- safe_numeric(
  evidence$FDR_siNRP1_5
)

evidence$significant_siNRP1_3 <- safe_logical(
  evidence$significant_siNRP1_3
)

evidence$significant_siNRP1_5 <- safe_logical(
  evidence$significant_siNRP1_5
)

evidence$in_signature_siNRP1_3 <- safe_logical(
  evidence$in_signature_siNRP1_3
)

evidence$in_signature_siNRP1_5 <- safe_logical(
  evidence$in_signature_siNRP1_5
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

evidence$modules <- safe_character(
  evidence$modules
)

cat("Input records:", nrow(evidence), "\n")
cat(
  "Unique genes:",
  length(unique(evidence$gene_symbol)),
  "\n"
)
cat(
  "Unique pathways:",
  length(unique(evidence$pathway_id)),
  "\n\n"
)


# ============================================================
# 4. REMOVE DUPLICATE PATHWAY-GENE RECORDS
# ============================================================

duplicate_count <- sum(
  duplicated(
    evidence[, c(
      "pathway_id",
      "gene_symbol"
    )]
  )
)

evidence <- evidence[
  !duplicated(
    evidence[, c(
      "pathway_id",
      "gene_symbol"
    )]
  ),
  ,
  drop = FALSE
]

cat(
  "Duplicate pathway-gene records removed:",
  duplicate_count,
  "\n"
)

cat(
  "Unique pathway-gene records:",
  nrow(evidence),
  "\n\n"
)


# ============================================================
# 5. DERIVE PATHWAY STATUS
# ============================================================

evidence$present_siNRP1_3 <-
  evidence$in_signature_siNRP1_3

evidence$present_siNRP1_5 <-
  evidence$in_signature_siNRP1_5

pathway_ids <- unique(
  evidence$pathway_id
)

pathway_records <- lapply(
  pathway_ids,
  function(current_id) {
    
    d <- evidence[
      evidence$pathway_id == current_id,
      ,
      drop = FALSE
    ]
    
    genes_3 <- unique(
      d$gene_symbol[
        d$present_siNRP1_3
      ]
    )
    
    genes_5 <- unique(
      d$gene_symbol[
        d$present_siNRP1_5
      ]
    )
    
    present_3 <- length(genes_3) > 0
    present_5 <- length(genes_5) > 0
    
    status <- if (
      present_3 && present_5
    ) {
      "SHARED"
    } else if (present_3) {
      "ONLY_siNRP1_3"
    } else if (present_5) {
      "ONLY_siNRP1_5"
    } else {
      "NOT_IN_SIGNATURE"
    }
    
    data.frame(
      pathway_id = current_id,
      term = d$term[1],
      ontology = d$ontology[1],
      group = d$group[1],
      pathway_status = status,
      signature_genes_siNRP1_3 = length(genes_3),
      signature_genes_siNRP1_5 = length(genes_5),
      stringsAsFactors = FALSE
    )
  }
)

pathway_summary <- do.call(
  rbind,
  pathway_records
)

rownames(pathway_summary) <- NULL


# ============================================================
# 6. GENE-LEVEL CONCORDANCE ANALYSIS
# ============================================================

pathway_statistics <- lapply(
  pathway_ids,
  function(current_id) {
    
    d <- evidence[
      evidence$pathway_id == current_id,
      ,
      drop = FALSE
    ]
    
    sig_3 <- d$significant_siNRP1_3
    sig_5 <- d$significant_siNRP1_5
    
    sig_both <- sig_3 & sig_5
    
    direction_3 <- d$direction_siNRP1_3
    direction_5 <- d$direction_siNRP1_5
    
    valid_direction <- sig_both &
      direction_3 %in% c("UP", "DOWN") &
      direction_5 %in% c("UP", "DOWN")
    
    concordant <- valid_direction &
      direction_3 == direction_5
    
    discordant <- valid_direction &
      direction_3 != direction_5
    
    modules <- parse_modules(
      d$modules
    )
    
    n_genes <- length(
      unique(d$gene_symbol)
    )
    
    n_sig_3 <- length(
      unique(d$gene_symbol[sig_3])
    )
    
    n_sig_5 <- length(
      unique(d$gene_symbol[sig_5])
    )
    
    n_both <- length(
      unique(d$gene_symbol[sig_both])
    )
    
    n_valid <- length(
      unique(d$gene_symbol[valid_direction])
    )
    
    n_concordant <- length(
      unique(d$gene_symbol[concordant])
    )
    
    n_discordant <- length(
      unique(d$gene_symbol[discordant])
    )
    
    module_gene_count <- length(
      unique(
        d$gene_symbol[
          nzchar(d$modules)
        ]
      )
    )
    
    data.frame(
      pathway_id = current_id,
      total_genes = n_genes,
      significant_genes_siNRP1_3 = n_sig_3,
      significant_genes_siNRP1_5 = n_sig_5,
      significant_both = n_both,
      evaluable_direction_genes = n_valid,
      concordant_genes = n_concordant,
      discordant_genes = n_discordant,
      concordance_fraction = safe_fraction(
        n_concordant,
        n_valid
      ),
      discordance_fraction = safe_fraction(
        n_discordant,
        n_valid
      ),
      module_associated_genes = module_gene_count,
      module_count = length(modules),
      modules = collapse_unique(modules),
      stringsAsFactors = FALSE
    )
  }
)

pathway_statistics <- do.call(
  rbind,
  pathway_statistics
)

rownames(pathway_statistics) <- NULL


# ============================================================
# 7. INTEGRATE PATHWAY-LEVEL RESULTS
# ============================================================

priority <- merge(
  pathway_summary,
  pathway_statistics,
  by = "pathway_id",
  all = TRUE,
  sort = FALSE
)

priority$significant_fraction_siNRP1_3 <-
  safe_fraction(
    priority$significant_genes_siNRP1_3,
    priority$total_genes
  )

priority$significant_fraction_siNRP1_5 <-
  safe_fraction(
    priority$significant_genes_siNRP1_5,
    priority$total_genes
  )

priority$module_gene_fraction <-
  safe_fraction(
    priority$module_associated_genes,
    priority$total_genes
  )

priority$directional_evidence_fraction <-
  safe_fraction(
    priority$evaluable_direction_genes,
    priority$total_genes
  )


# ============================================================
# 8. PATHWAY PRIORITY SCORE
# ============================================================

# Transparent weighted score:
#
# 30% = gene support
# 25% = significant gene fraction
# 25% = directional concordance
# 20% = module-associated gene fraction
#
# Each component is normalized across the pathways.
#
# This score is a descriptive ranking metric.
# It is not a statistical significance measure.

priority$score_gene_support <- minmax_score(
  log1p(priority$total_genes)
)

priority$score_significant_support <- minmax_score(
  (
    priority$significant_genes_siNRP1_3 +
      priority$significant_genes_siNRP1_5
  ) /
    2
)

priority$score_concordance <- minmax_score(
  priority$concordance_fraction
)

priority$score_module_support <- minmax_score(
  priority$module_gene_fraction
)

priority$priority_score <-
  0.30 * priority$score_gene_support +
  0.25 * priority$score_significant_support +
  0.25 * priority$score_concordance +
  0.20 * priority$score_module_support

priority$priority_class <- ifelse(
  priority$priority_score >= 75,
  "HIGH_PRIORITY",
  ifelse(
    priority$priority_score >= 50,
    "MEDIUM_PRIORITY",
    "EXPLORATORY"
  )
)

priority$priority_score <- round(
  priority$priority_score,
  2
)

priority$concordance_fraction <- round(
  priority$concordance_fraction,
  4
)

priority$discordance_fraction <- round(
  priority$discordance_fraction,
  4
)

priority$significant_fraction_siNRP1_3 <- round(
  priority$significant_fraction_siNRP1_3,
  4
)

priority$significant_fraction_siNRP1_5 <- round(
  priority$significant_fraction_siNRP1_5,
  4
)

priority$module_gene_fraction <- round(
  priority$module_gene_fraction,
  4
)

priority$directional_evidence_fraction <- round(
  priority$directional_evidence_fraction,
  4
)


# ============================================================
# 9. KNOCKDOWN-SPECIFIC EVIDENCE
# ============================================================

priority$stronger_gene_support <- ifelse(
  priority$significant_genes_siNRP1_3 >
    priority$significant_genes_siNRP1_5,
  "siNRP1_3",
  ifelse(
    priority$significant_genes_siNRP1_5 >
      priority$significant_genes_siNRP1_3,
    "siNRP1_5",
    "BALANCED"
  )
)

priority$evidence_difference <-
  priority$significant_genes_siNRP1_5 -
  priority$significant_genes_siNRP1_3

priority <- priority[
  order(
    -priority$priority_score,
    -priority$significant_both,
    -priority$total_genes
  ),
  ,
  drop = FALSE
]

priority$priority_rank <- seq_len(
  nrow(priority)
)

rownames(priority) <- NULL


# ============================================================
# 10. SHARED PATHWAY ANALYSIS
# ============================================================

shared_priority <- priority[
  priority$pathway_status == "SHARED",
  ,
  drop = FALSE
]

shared_priority <- shared_priority[
  order(
    -shared_priority$priority_score,
    -shared_priority$significant_both
  ),
  ,
  drop = FALSE
]

shared_priority$shared_rank <- seq_len(
  nrow(shared_priority)
)

rownames(shared_priority) <- NULL


# ============================================================
# 11. GENE-LEVEL DISCORDANCE TABLE
# ============================================================

discordance_records <- evidence[
  evidence$significant_both &
    evidence$direction_siNRP1_3 %in% c("UP", "DOWN") &
    evidence$direction_siNRP1_5 %in% c("UP", "DOWN"),
  ,
  drop = FALSE
]

discordance_records$concordance_status <- ifelse(
  discordance_records$direction_siNRP1_3 ==
    discordance_records$direction_siNRP1_5,
  "CONCORDANT",
  "DISCORDANT"
)

discordance_records$logFC_difference <-
  discordance_records$logFC_siNRP1_5 -
  discordance_records$logFC_siNRP1_3

discordance_records$absolute_logFC_difference <- abs(
  discordance_records$logFC_difference
)

discordance_records <- merge(
  discordance_records,
  priority[, c(
    "pathway_id",
    "pathway_status",
    "priority_score",
    "priority_class"
  )],
  by = "pathway_id",
  all.x = TRUE,
  sort = FALSE
)

discordance_records <- discordance_records[
  order(
    discordance_records$concordance_status,
    -discordance_records$absolute_logFC_difference
  ),
  ,
  drop = FALSE
]

discordance_output <- discordance_records[, c(
  "gene_symbol",
  "pathway_id",
  "term",
  "ontology",
  "group",
  "pathway_status",
  "direction_siNRP1_3",
  "direction_siNRP1_5",
  "logFC_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_3",
  "FDR_siNRP1_5",
  "logFC_difference",
  "absolute_logFC_difference",
  "concordance_status",
  "modules",
  "priority_score",
  "priority_class"
)]


# ============================================================
# 12. MODULE SUPPORT SUMMARY
# ============================================================

module_records <- list()

for (i in seq_len(nrow(evidence))) {
  
  current_modules <- parse_modules(
    evidence$modules[i]
  )
  
  if (length(current_modules) == 0) {
    next
  }
  
  for (current_module in current_modules) {
    
    module_records[[length(module_records) + 1]] <-
      data.frame(
        pathway_id = evidence$pathway_id[i],
        gene_symbol = evidence$gene_symbol[i],
        module = current_module,
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
  
  module_data <- do.call(
    rbind,
    module_records
  )
  
  module_data <- module_data[
    !duplicated(
      module_data[, c(
        "pathway_id",
        "gene_symbol",
        "module"
      )]
    ),
    ,
    drop = FALSE
  ]
  
  module_groups <- split(
    module_data,
    module_data$module
  )
  
  module_summary <- lapply(
    names(module_groups),
    function(module_name) {
      
      d <- module_groups[[module_name]]
      
      sig_both <- d$significant_both
      
      valid_direction <- sig_both &
        d$direction_siNRP1_3 %in% c("UP", "DOWN") &
        d$direction_siNRP1_5 %in% c("UP", "DOWN")
      
      concordant <- valid_direction &
        d$direction_siNRP1_3 ==
        d$direction_siNRP1_5
      
      discordant <- valid_direction &
        d$direction_siNRP1_3 !=
        d$direction_siNRP1_5
      
      data.frame(
        module = module_name,
        pathway_count = length(
          unique(d$pathway_id)
        ),
        unique_gene_count = length(
          unique(d$gene_symbol)
        ),
        significant_siNRP1_3 = length(
          unique(
            d$gene_symbol[
              d$significant_siNRP1_3
            ]
          )
        ),
        significant_siNRP1_5 = length(
          unique(
            d$gene_symbol[
              d$significant_siNRP1_5
            ]
          )
        ),
        significant_both = length(
          unique(
            d$gene_symbol[
              sig_both
            ]
          )
        ),
        concordant_genes = length(
          unique(
            d$gene_symbol[
              concordant
            ]
          )
        ),
        discordant_genes = length(
          unique(
            d$gene_symbol[
              discordant
            ]
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  module_summary <- do.call(
    rbind,
    module_summary
  )
  
  module_summary$directional_concordance <- round(
    safe_fraction(
      module_summary$concordant_genes,
      module_summary$concordant_genes +
        module_summary$discordant_genes
    ),
    4
  )
  
  module_summary <- module_summary[
    order(
      -module_summary$significant_both,
      -module_summary$unique_gene_count
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  module_summary <- data.frame(
    module = character(0),
    pathway_count = integer(0),
    unique_gene_count = integer(0),
    significant_siNRP1_3 = integer(0),
    significant_siNRP1_5 = integer(0),
    significant_both = integer(0),
    concordant_genes = integer(0),
    discordant_genes = integer(0),
    directional_concordance = numeric(0),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 13. PATHWAY STATUS SUMMARY
# ============================================================

status_summary <- as.data.frame(
  table(priority$pathway_status),
  stringsAsFactors = FALSE
)

colnames(status_summary) <- c(
  "pathway_status",
  "pathway_count"
)

status_summary <- status_summary[
  order(-status_summary$pathway_count),
  ,
  drop = FALSE
]


# ============================================================
# 14. PRIORITY CLASS SUMMARY
# ============================================================

priority_class_summary <- as.data.frame(
  table(priority$priority_class),
  stringsAsFactors = FALSE
)

colnames(priority_class_summary) <- c(
  "priority_class",
  "pathway_count"
)


# ============================================================
# 15. AUDIT
# ============================================================

audit <- data.frame(
  metric = c(
    "input_file",
    "input_records",
    "records_after_duplicate_removal",
    "duplicate_records_removed",
    "unique_genes",
    "unique_pathways",
    "shared_pathways",
    "only_siNRP1_3_pathways",
    "only_siNRP1_5_pathways",
    "pathway_gene_records_significant_both",
    "concordant_pathway_gene_records",
    "discordant_pathway_gene_records",
    "pathways_high_priority",
    "pathways_medium_priority",
    "pathways_exploratory",
    "module_categories"
  ),
  value = c(
    INPUT_FILE,
    nrow(read.csv(
      INPUT_FILE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )),
    nrow(evidence),
    duplicate_count,
    length(unique(evidence$gene_symbol)),
    length(unique(evidence$pathway_id)),
    sum(priority$pathway_status == "SHARED"),
    sum(priority$pathway_status == "ONLY_siNRP1_3"),
    sum(priority$pathway_status == "ONLY_siNRP1_5"),
    sum(evidence$significant_both),
    sum(
      discordance_records$concordance_status ==
        "CONCORDANT"
    ),
    sum(
      discordance_records$concordance_status ==
        "DISCORDANT"
    ),
    sum(priority$priority_class == "HIGH_PRIORITY"),
    sum(priority$priority_class == "MEDIUM_PRIORITY"),
    sum(priority$priority_class == "EXPLORATORY"),
    nrow(module_summary)
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 16. EXPORT TABLES
# ============================================================

write_csv(
  priority,
  file.path(
    TABLE_DIR,
    "pathway_priority_ranking.csv"
  )
)

write_csv(
  shared_priority,
  file.path(
    TABLE_DIR,
    "shared_pathway_priority.csv"
  )
)

write_csv(
  discordance_output,
  file.path(
    TABLE_DIR,
    "pathway_gene_concordance_discordance.csv"
  )
)

write_csv(
  module_summary,
  file.path(
    TABLE_DIR,
    "module_support_summary.csv"
  )
)

write_csv(
  status_summary,
  file.path(
    TABLE_DIR,
    "pathway_status_summary.csv"
  )
)

write_csv(
  priority_class_summary,
  file.path(
    TABLE_DIR,
    "priority_class_summary.csv"
  )
)

write_csv(
  audit,
  file.path(
    AUDIT_DIR,
    "script36_audit.csv"
  )
)


# ============================================================
# 17. REPORT
# ============================================================

report_file <- file.path(
  REPORT_DIR,
  "shared_pathway_priority_report.txt"
)

sink(report_file)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 36 - SHARED PATHWAY PRIORITY ANALYSIS\n")
cat("============================================================\n\n")

cat("INPUT\n")
cat("Pathway-gene evidence records:", nrow(evidence), "\n")
cat("Unique genes:", length(unique(evidence$gene_symbol)), "\n")
cat("Unique pathways:", length(unique(evidence$pathway_id)), "\n\n")

cat("PATHWAY STATUS\n")
print(status_summary, row.names = FALSE)

cat("\nPRIORITY CLASS DISTRIBUTION\n")
print(priority_class_summary, row.names = FALSE)

cat("\nSHARED PATHWAYS\n")
cat("Number of shared pathways:", nrow(shared_priority), "\n\n")

if (nrow(shared_priority) > 0) {
  
  print(
    shared_priority[, c(
      "shared_rank",
      "term",
      "priority_score",
      "priority_class",
      "total_genes",
      "significant_both",
      "concordant_genes",
      "discordant_genes",
      "modules"
    )],
    row.names = FALSE
  )
}

cat("\nGENE-LEVEL DIRECTIONAL EVIDENCE\n")

cat(
  "Significant in both contrasts:",
  sum(evidence$significant_both),
  "\n"
)

cat(
  "Concordant records:",
  sum(
    discordance_records$concordance_status ==
      "CONCORDANT"
  ),
  "\n"
)

cat(
  "Discordant records:",
  sum(
    discordance_records$concordance_status ==
      "DISCORDANT"
  ),
  "\n"
)

cat("\nMODULE SUPPORT\n")

if (nrow(module_summary) > 0) {
  print(
    module_summary,
    row.names = FALSE
  )
} else {
  cat("No module annotations available.\n")
}

cat("\nTOP 15 PRIORITIZED PATHWAYS\n")

top_pathways <- head(
  priority,
  15
)

print(
  top_pathways[, c(
    "priority_rank",
    "term",
    "pathway_status",
    "priority_score",
    "priority_class",
    "significant_genes_siNRP1_3",
    "significant_genes_siNRP1_5",
    "concordant_genes",
    "discordant_genes"
  )],
  row.names = FALSE
)

cat("\nMETHODOLOGICAL NOTES\n")
cat(
  "- Pathway status reflects signature membership in Script 35.\n"
)
cat(
  "- Gene counts are unique within each pathway.\n"
)
cat(
  "- A gene may contribute to multiple pathways.\n"
)
cat(
  "- Concordance is evaluated only for genes significant in both contrasts with valid directions.\n"
)
cat(
  "- The priority score is descriptive and uses the documented weighted components.\n"
)
cat(
  "- Pathway prioritization does not demonstrate pathway activation or inhibition.\n"
)
cat(
  "- No new enrichment tests or differential expression analyses were performed.\n"
)

cat("\nOUTPUT DIRECTORY\n")
cat(OUTPUT_DIR, "\n")

cat("\n============================================================\n")
cat("END OF SCRIPT 36\n")
cat("============================================================\n")

sink()


# ============================================================
# 18. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 36 COMPLETED\n")
cat("============================================================\n\n")

cat("Pathways analyzed:", nrow(priority), "\n")
cat("Shared pathways:", nrow(shared_priority), "\n")

cat(
  "High priority:",
  sum(priority$priority_class == "HIGH_PRIORITY"),
  "\n"
)

cat(
  "Medium priority:",
  sum(priority$priority_class == "MEDIUM_PRIORITY"),
  "\n"
)

cat(
  "Exploratory:",
  sum(priority$priority_class == "EXPLORATORY"),
  "\n"
)

cat("\nOUTPUT FILES\n")

cat(
  "1. pathway_priority_ranking.csv\n",
  "2. shared_pathway_priority.csv\n",
  "3. pathway_gene_concordance_discordance.csv\n",
  "4. module_support_summary.csv\n",
  "5. pathway_status_summary.csv\n",
  "6. priority_class_summary.csv\n",
  "7. script36_audit.csv\n",
  "8. shared_pathway_priority_report.txt\n",
  sep = ""
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n============================================================\n")
