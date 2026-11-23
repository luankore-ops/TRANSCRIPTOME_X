
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 51 - FINAL BIOLOGICAL CONCLUSIONS AND RESEARCH HYPOTHESES
#
# Objective:
# Integrate final candidate prioritization into a concise,
# scientifically cautious biological conclusions report.
#
# Input:
# results/final_candidate_biological_evidence_assessment/
#
# Output:
# results/final_biological_conclusions/
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 51 - FINAL BIOLOGICAL CONCLUSIONS AND HYPOTHESES\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

project_dir <- paste0(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/",
  "TRANSCRIPTOME_X"
)

input_dir <- file.path(
  project_dir,
  "results",
  "final_candidate_biological_evidence_assessment"
)

input_tables <- file.path(input_dir, "tables")

output_dir <- file.path(
  project_dir,
  "results",
  "final_biological_conclusions"
)

tables_dir <- file.path(output_dir, "tables")
reports_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (folder in c(tables_dir, reports_dir, audit_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

required_files <- c(
  "final_candidate_biological_assessment.csv",
  "final_candidate_pathway_landscape.csv",
  "final_candidate_module_landscape.csv",
  "pathway_module_candidate_convergence.csv",
  "final_candidate_evidence_classification.csv",
  "final_candidate_functional_support.csv"
)

missing_files <- required_files[
  !file.exists(file.path(input_tables, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing Script 50 outputs:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------

read_input <- function(filename) {
  read.csv(
    file.path(input_tables, filename),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

assessment <- read_input(
  "final_candidate_biological_assessment.csv"
)

pathway_landscape <- read_input(
  "final_candidate_pathway_landscape.csv"
)

module_landscape <- read_input(
  "final_candidate_module_landscape.csv"
)

convergence <- read_input(
  "pathway_module_candidate_convergence.csv"
)

classification <- read_input(
  "final_candidate_evidence_classification.csv"
)

functional_support <- read_input(
  "final_candidate_functional_support.csv"
)

cat("Input datasets loaded.\n")
cat("Candidates:", nrow(assessment), "\n")

# ------------------------------------------------------------
# 3. LOAD ASSOCIATION TABLES FROM SCRIPT 49
# ------------------------------------------------------------

script49_dir <- file.path(
  project_dir,
  "results",
  "integrated_biological_evidence_synthesis",
  "tables"
)

pathway_file <- file.path(
  script49_dir,
  "candidate_pathway_evidence.csv"
)

module_file <- file.path(
  script49_dir,
  "candidate_module_evidence.csv"
)

if (!file.exists(pathway_file) || !file.exists(module_file)) {
  stop("Script 49 association tables not found.")
}

pathway_edges <- read.csv(
  pathway_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

module_edges <- read.csv(
  module_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 4. HELPER FUNCTIONS
# ------------------------------------------------------------

find_column <- function(data, candidates, required = FALSE) {
  
  available <- names(data)
  
  indices <- match(
    tolower(candidates),
    tolower(available)
  )
  
  indices <- indices[!is.na(indices)]
  
  if (length(indices) > 0) {
    return(available[indices[1]])
  }
  
  if (required) {
    stop(
      "Required column not found: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns: ",
      paste(available, collapse = ", ")
    )
  }
  
  NULL
}

unique_values <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      nzchar(trimws(x))
  ]
  
  unique(x)
}

collapse_values <- function(x, limit = 5) {
  
  values <- unique_values(x)
  
  if (length(values) == 0) {
    return("Not available")
  }
  
  shown <- head(values, limit)
  
  result <- paste(shown, collapse = "; ")
  
  if (length(values) > limit) {
    result <- paste0(
      result,
      "; and ",
      length(values) - limit,
      " additional associations"
    )
  }
  
  result
}

safe_write <- function(data, filename) {
  
  write.csv(
    data,
    file.path(tables_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

# ------------------------------------------------------------
# 5. IDENTIFY ASSOCIATION COLUMNS
# ------------------------------------------------------------

gene_assessment_col <- find_column(
  assessment,
  c("gene_symbol"),
  required = TRUE
)

gene_pathway_col <- find_column(
  pathway_edges,
  c("gene_symbol", "candidate_gene", "gene", "symbol"),
  required = TRUE
)

gene_module_col <- find_column(
  module_edges,
  c("gene_symbol", "candidate_gene", "gene", "symbol"),
  required = TRUE
)

pathway_col <- find_column(
  pathway_edges,
  c("pathway_name", "pathway", "term_name", "pathway_id", "term"),
  required = TRUE
)

module_col <- find_column(
  module_edges,
  c("module_name", "module", "gene_module"),
  required = TRUE
)

pathway_edges$gene_symbol <- as.character(
  pathway_edges[[gene_pathway_col]]
)

pathway_edges$pathway_label <- as.character(
  pathway_edges[[pathway_col]]
)

module_edges$gene_symbol <- as.character(
  module_edges[[gene_module_col]]
)

module_edges$module_label <- as.character(
  module_edges[[module_col]]
)

assessment$gene_symbol <- as.character(
  assessment[[gene_assessment_col]]
)

# ------------------------------------------------------------
# 6. DEFINE PRIORITY CANDIDATES
# ------------------------------------------------------------

rank_col <- find_column(
  assessment,
  c("final_assessment_rank"),
  required = TRUE
)

assessment[[rank_col]] <- suppressWarnings(
  as.numeric(assessment[[rank_col]])
)

assessment <- assessment[
  order(assessment[[rank_col]]),
  ,
  drop = FALSE
]

top_candidates <- head(assessment, 20)

# Keep only associations belonging to assessed candidates.
pathway_edges <- pathway_edges[
  pathway_edges$gene_symbol %in% assessment$gene_symbol,
  ,
  drop = FALSE
]

module_edges <- module_edges[
  module_edges$gene_symbol %in% assessment$gene_symbol,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 7. CANDIDATE-LEVEL FUNCTIONAL EVIDENCE
# ------------------------------------------------------------

candidate_evidence <- data.frame(
  final_assessment_rank = integer(),
  gene_symbol = character(),
  pathway_count = integer(),
  module_count = integer(),
  associated_pathways = character(),
  associated_modules = character(),
  evidence_classification = character(),
  interpretation = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(top_candidates))) {
  
  gene <- top_candidates$gene_symbol[i]
  
  gene_pathways <- pathway_edges[
    pathway_edges$gene_symbol == gene,
    ,
    drop = FALSE
  ]
  
  gene_modules <- module_edges[
    module_edges$gene_symbol == gene,
    ,
    drop = FALSE
  ]
  
  pathways <- unique_values(gene_pathways$pathway_label)
  modules <- unique_values(gene_modules$module_label)
  
  classification_col <- find_column(
    top_candidates,
    c("evidence_assessment")
  )
  
  classification_value <- if (!is.null(classification_col)) {
    as.character(top_candidates[[classification_col]][i])
  } else {
    "Not available"
  }
  
  interpretation <- if (
    length(pathways) > 0 &&
    length(modules) > 0
  ) {
    paste(
      "The candidate is represented in both pathway",
      "and module associations. These observations support",
      "functional connectivity within the analyzed dataset."
    )
  } else if (length(pathways) > 0) {
    paste(
      "The candidate has pathway associations in the",
      "integrated evidence dataset."
    )
  } else if (length(modules) > 0) {
    paste(
      "The candidate has module associations in the",
      "integrated evidence dataset."
    )
  } else {
    "No pathway or module associations were identified."
  }
  
  candidate_evidence <- rbind(
    candidate_evidence,
    data.frame(
      final_assessment_rank =
        top_candidates[[rank_col]][i],
      gene_symbol = gene,
      pathway_count = length(pathways),
      module_count = length(modules),
      associated_pathways = collapse_values(pathways, 8),
      associated_modules = collapse_values(modules, 8),
      evidence_classification = classification_value,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------
# 8. NRP1 FOCUSED EVIDENCE SUMMARY
# ------------------------------------------------------------

nrp1_present <- "NRP1" %in% assessment$gene_symbol

if (nrp1_present) {
  
  nrp1_assessment <- assessment[
    assessment$gene_symbol == "NRP1",
    ,
    drop = FALSE
  ]
  
  nrp1_pathways <- unique_values(
    pathway_edges$pathway_label[
      pathway_edges$gene_symbol == "NRP1"
    ]
  )
  
  nrp1_modules <- unique_values(
    module_edges$module_label[
      module_edges$gene_symbol == "NRP1"
    ]
  )
  
  nrp1_summary <- data.frame(
    gene_symbol = "NRP1",
    final_assessment_rank = nrp1_assessment[[rank_col]][1],
    pathway_count = length(nrp1_pathways),
    module_count = length(nrp1_modules),
    associated_pathways = collapse_values(nrp1_pathways, 20),
    associated_modules = collapse_values(nrp1_modules, 10),
    stringsAsFactors = FALSE
  )
  
} else {
  
  nrp1_summary <- data.frame(
    gene_symbol = "NRP1",
    final_assessment_rank = NA_real_,
    pathway_count = 0,
    module_count = 0,
    associated_pathways = "NRP1 not present in candidate assessment",
    associated_modules = "NRP1 not present in candidate assessment",
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 9. RESEARCH HYPOTHESES
# ------------------------------------------------------------

hypotheses <- data.frame(
  hypothesis_id = character(),
  candidate_gene = character(),
  evidence_basis = character(),
  research_hypothesis = character(),
  validation_strategy = character(),
  evidence_status = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(candidate_evidence))) {
  
  gene <- candidate_evidence$gene_symbol[i]
  
  pathway_text <- candidate_evidence$associated_pathways[i]
  module_text <- candidate_evidence$associated_modules[i]
  
  hypothesis_text <- paste0(
    gene,
    " may participate in biological processes represented by ",
    "its pathway and module associations in the NRP1 knockdown ",
    "transcriptomic dataset. Its functional relevance should be ",
    "tested experimentally."
  )
  
  validation_text <- paste(
    "Independent expression validation, followed by targeted",
    "perturbation and functional assays relevant to the",
    "associated biological processes."
  )
  
  hypotheses <- rbind(
    hypotheses,
    data.frame(
      hypothesis_id = paste0("H", i),
      candidate_gene = gene,
      evidence_basis = paste0(
        "Pathways: ", pathway_text,
        " | Modules: ", module_text
      ),
      research_hypothesis = hypothesis_text,
      validation_strategy = validation_text,
      evidence_status = "COMPUTATIONAL HYPOTHESIS",
      stringsAsFactors = FALSE
    )
  )
}

# Dedicated NRP1 research hypothesis.
if (nrp1_present) {
  
  nrp1_hypothesis <- data.frame(
    hypothesis_id = "H_NRP1",
    candidate_gene = "NRP1",
    evidence_basis = paste0(
      "NRP1 ranked ",
      nrp1_summary$final_assessment_rank,
      " in the final assessment, with ",
      nrp1_summary$pathway_count,
      " pathway associations and ",
      nrp1_summary$module_count,
      " module associations."
    ),
    research_hypothesis = paste(
      "NRP1 knockdown may be associated with transcriptional",
      "changes involving multiple functionally annotated",
      "biological processes. NRP1 may therefore represent",
      "a relevant candidate for investigating the molecular",
      "response observed in this cellular model."
    ),
    validation_strategy = paste(
      "Validate NRP1 knockdown efficiency, replicate the",
      "transcriptional findings in independent experiments,",
      "and test selected downstream candidates using targeted",
      "molecular and cellular assays."
    ),
    evidence_status = "COMPUTATIONAL HYPOTHESIS - PRIORITY CANDIDATE",
    stringsAsFactors = FALSE
  )
  
  hypotheses <- rbind(
    nrp1_hypothesis,
    hypotheses
  )
}

# ------------------------------------------------------------
# 10. FINAL BIOLOGICAL CONCLUSIONS
# ------------------------------------------------------------

top_gene_names <- paste(
  head(assessment$gene_symbol, 10),
  collapse = ", "
)

top_candidate_pathways <- head(
  pathway_landscape,
  10
)

top_candidate_modules <- head(
  module_landscape,
  10
)

nrp1_rank <- if (nrp1_present) {
  nrp1_summary$final_assessment_rank
} else {
  NA
}

nrp1_pathway_count <- if (nrp1_present) {
  nrp1_summary$pathway_count
} else {
  0
}

nrp1_module_count <- if (nrp1_present) {
  nrp1_summary$module_count
} else {
  0
}

conclusion_text <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "FINAL BIOLOGICAL CONCLUSIONS AND RESEARCH HYPOTHESES",
  "============================================================",
  "",
  paste("Report generated:", Sys.time()),
  "",
  "1. STUDY CONTEXT",
  "------------------------------------------------------------",
  paste(
    "This project investigated transcriptomic changes associated",
    "with NRP1 knockdown in the MDA-MB-231 claudin-low breast",
    "cancer cell model, using the analyzed GEO dataset GSE266566."
  ),
  "",
  "2. INTEGRATED CANDIDATE PRIORITIZATION",
  "------------------------------------------------------------",
  paste("Candidates assessed:", nrow(assessment)),
  paste("Top ten candidates:", top_gene_names),
  paste(
    "Distinct pathway annotations represented:",
    nrow(pathway_landscape)
  ),
  paste(
    "Distinct modules represented:",
    nrow(module_landscape)
  ),
  "",
  "The integrated computational analysis prioritized a set of",
  "candidate genes supported by previously defined evidence",
  "layers, including differential expression, network context,",
  "pathway associations and module associations.",
  "",
  "3. PRIORITY CANDIDATES",
  "------------------------------------------------------------",
  capture.output(
    print(
      candidate_evidence[
        ,
        c(
          "final_assessment_rank",
          "gene_symbol",
          "pathway_count",
          "module_count",
          "evidence_classification"
        )
      ],
      row.names = FALSE
    )
  ),
  "",
  "The candidates should be interpreted as computationally",
  "prioritized genes rather than experimentally confirmed",
  "regulators or causal mediators.",
  "",
  "4. NRP1-FOCUSED INTERPRETATION",
  "------------------------------------------------------------",
  if (nrp1_present) {
    paste(
      "NRP1 was included in the final candidate assessment",
      "at rank", nrp1_rank, "with",
      nrp1_pathway_count, "pathway associations and",
      nrp1_module_count, "module associations."
    )
  } else {
    "NRP1 was not found in the final candidate assessment."
  },
  "",
  paste(
    "These associations provide a basis for investigating",
    "the broader transcriptional response accompanying",
    "NRP1 knockdown. They do not establish that NRP1 directly",
    "regulates each associated pathway or module."
  ),
  "",
  "5. BIOLOGICAL INTERPRETATION",
  "------------------------------------------------------------",
  paste(
    "The candidate landscape suggests that the observed",
    "transcriptional response involves genes connected to",
    "multiple annotated biological processes."
  ),
  paste(
    "Genes such as THBS1, CCN3, NRP1 and TGFBR3 are among",
    "the highest-ranked candidates in the integrated assessment."
  ),
  paste(
    "Their prioritization provides a rational basis for",
    "selecting targets for subsequent experimental investigation."
  ),
  "",
  "6. RESEARCH HYPOTHESES",
  "------------------------------------------------------------",
  paste(
    "H1: NRP1 knockdown is associated with a broader",
    "transcriptional response involving multiple functional",
    "gene associations in this cellular model."
  ),
  paste(
    "H2: A subset of prioritized genes, including THBS1, CCN3",
    "and TGFBR3, may represent relevant components of this",
    "transcriptional response and warrant independent validation."
  ),
  paste(
    "H3: Functional relationships among prioritized genes",
    "may help identify biological processes for targeted",
    "follow-up experiments."
  ),
  "",
  "These hypotheses are generated from computational evidence",
  "and require independent experimental testing.",
  "",
  "7. RECOMMENDED EXPERIMENTAL VALIDATION",
  "------------------------------------------------------------",
  "1. Confirm NRP1 knockdown efficiency at RNA and protein level.",
  "2. Validate selected candidate genes using independent qPCR.",
  "3. Confirm selected protein-level changes where appropriate.",
  "4. Perform targeted functional assays based on the strongest",
  "   pathway and module associations.",
  "5. Use independent biological replicates and, where possible,",
  "   additional cellular models.",
  "",
  "8. LIMITATIONS",
  "------------------------------------------------------------",
  "- The conclusions are based on a computational analysis of",
  "  a specific cellular model and dataset.",
  "- Pathway associations do not demonstrate pathway activation.",
  "- Network relationships do not establish direct regulation.",
  "- Candidate ranking is not equivalent to experimental validation.",
  "- Evidence layers partly reuse indicators generated in previous",
  "  integrated analyses and are not fully independent.",
  "- The number of pathway associations is not a measure of effect",
  "  size or biological importance.",
  "- Findings should not be generalized directly to patients or",
  "  other breast cancer models without further validation.",
  "",
  "9. OVERALL CONCLUSION",
  "------------------------------------------------------------",
  paste(
    "TRANSCRIPTOME-X provides an integrated computational",
    "framework for prioritizing genes and functional associations",
    "related to the transcriptional response observed following",
    "NRP1 knockdown in MDA-MB-231 cells."
  ),
  paste(
    "The resulting candidate panel, particularly the highest-ranked",
    "genes, provides a focused starting point for experimental",
    "validation and future mechanistic studies."
  ),
  "",
  "The results support candidate prioritization and hypothesis",
  "generation, not definitive causal conclusions.",
  "",
  "============================================================",
  "END OF REPORT",
  "============================================================"
)

writeLines(
  conclusion_text,
  file.path(
    reports_dir,
    "TRANSCRIPTOME_X_final_biological_conclusions.txt"
  )
)

# ------------------------------------------------------------
# 11. EXPORT TABLES
# ------------------------------------------------------------

safe_write(
  candidate_evidence,
  "final_priority_candidate_biological_evidence.csv"
)

safe_write(
  nrp1_summary,
  "nrp1_focused_evidence_summary.csv"
)

safe_write(
  hypotheses,
  "final_research_hypotheses.csv"
)

safe_write(
  top_candidates,
  "top_20_final_candidate_assessment.csv"
)

safe_write(
  pathway_landscape,
  "final_pathway_landscape.csv"
)

safe_write(
  module_landscape,
  "final_module_landscape.csv"
)

safe_write(
  convergence,
  "final_pathway_module_convergence.csv"
)

# ------------------------------------------------------------
# 12. CONSISTENCY CHECKS
# ------------------------------------------------------------

checks <- data.frame(
  check = c(
    "Candidate identifiers are unique",
    "All top candidates belong to assessment",
    "All candidate evidence rows have gene identifiers",
    "NRP1 summary generated",
    "Research hypotheses generated",
    "All hypothesis records have candidate identifiers",
    "Final report generated"
  ),
  passed = c(
    !anyDuplicated(assessment$gene_symbol),
    all(top_candidates$gene_symbol %in% assessment$gene_symbol),
    all(!is.na(candidate_evidence$gene_symbol)),
    nrow(nrp1_summary) == 1,
    nrow(hypotheses) > 0,
    all(!is.na(hypotheses$candidate_gene)),
    file.exists(
      file.path(
        reports_dir,
        "TRANSCRIPTOME_X_final_biological_conclusions.txt"
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  checks,
  file.path(audit_dir, "script51_consistency_checks.csv"),
  row.names = FALSE
)

if (!all(checks$passed)) {
  print(checks)
  stop("Script 51 consistency checks failed.")
}

# ------------------------------------------------------------
# 13. AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "assessed_candidates",
    "top_candidates",
    "pathway_associations",
    "module_associations",
    "research_hypotheses",
    "consistency_checks_passed"
  ),
  value = c(
    nrow(assessment),
    nrow(top_candidates),
    nrow(pathway_edges),
    nrow(module_edges),
    nrow(hypotheses),
    sum(checks$passed)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(audit_dir, "script51_audit.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "SCRIPT 51 EXECUTION LOG",
    paste("Execution:", Sys.time()),
    paste("Candidates:", nrow(assessment)),
    paste("Hypotheses:", nrow(hypotheses)),
    paste("Checks passed:", sum(checks$passed), "/", nrow(checks)),
    "Status: SUCCESS"
  ),
  file.path(audit_dir, "script51.log")
)

# ------------------------------------------------------------
# 14. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 51 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Candidates assessed:", nrow(assessment), "\n")
cat("Top candidates:", nrow(top_candidates), "\n")
cat("Research hypotheses:", nrow(hypotheses), "\n")
cat("Consistency checks:", sum(checks$passed), "/", nrow(checks), "\n\n")

cat("NRP1 FOCUSED SUMMARY\n")
print(nrp1_summary, row.names = FALSE)

cat("\nTOP PRIORITY CANDIDATES\n")
print(
  candidate_evidence[
    ,
    c(
      "final_assessment_rank",
      "gene_symbol",
      "pathway_count",
      "module_count"
    )
  ],
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\n============================================================\n")
