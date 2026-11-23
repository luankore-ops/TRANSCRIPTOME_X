# ============================================================
# TRANSCRIPTOME-X
# Script 34: Integrated pathway-module-DE analysis
#
# Integra i risultati validati dello script 33.
# Non esegue nuovi test DE o pathway enrichment.
# Le sovrapposizioni pathway-modulo sono descrittive.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

find_project_root <- function() {
  candidates <- c(getwd())

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    candidates <- c(candidates, dirname(dirname(normalizePath(
      script_path, winslash = "/", mustWork = FALSE
    ))))
  }

  for (candidate in unique(candidates)) {
    required <- file.path(
      candidate,
      "results",
      "integrated_gene_module_analysis",
      "tables",
      "module_integrated_gene_comparison.csv"
    )

    if (file.exists(required)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Cartella TRANSCRIPTOME_X non individuata. ",
    "Eseguire Rscript dalla directory principale del progetto."
  )
}

PROJECT_ROOT <- find_project_root()

INPUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "integrated_gene_module_analysis"
)

OUTPUT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  "integrated_pathway_module_analysis"
)

OUTPUT_TABLES <- file.path(OUTPUT_ROOT, "tables")
OUTPUT_AUDIT <- file.path(OUTPUT_ROOT, "audit")
OUTPUT_REPORTS <- file.path(OUTPUT_ROOT, "reports")

for (directory in c(
  OUTPUT_ROOT, OUTPUT_TABLES, OUTPUT_AUDIT, OUTPUT_REPORTS
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# 2. INPUT FILES
# ------------------------------------------------------------

INPUT_FILES <- list(
  integrated = file.path(
    INPUT_ROOT, "tables", "module_integrated_gene_comparison.csv"
  ),
  overlap = file.path(
    INPUT_ROOT, "tables", "signature_module_gene_overlap.csv"
  ),
  definitions = file.path(
    INPUT_ROOT, "tables", "module_gene_definitions.csv"
  )
)

missing_files <- names(INPUT_FILES)[
  !file.exists(unlist(INPUT_FILES, use.names = FALSE))
]

if (length(missing_files) > 0) {
  stop(
    "Input mancanti: ",
    paste(missing_files, collapse = ", ")
  )
}

read_input <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

integrated <- read_input(INPUT_FILES$integrated)
overlap <- read_input(INPUT_FILES$overlap)
definitions <- read_input(INPUT_FILES$definitions)
# ------------------------------------------------------------
# 3. REQUIRED COLUMNS AND BASIC NORMALIZATION
# ------------------------------------------------------------

required_columns <- list(
  integrated = c(
    "gene_symbol", "module",
    "logFC_siNRP1_3", "FDR_siNRP1_3",
    "logFC_siNRP1_5", "FDR_siNRP1_5",
    "significant_siNRP1_3", "significant_siNRP1_5",
    "direction_concordance", "joint_class"
  ),
  overlap = c(
    "gene_symbol", "signature_row", "contrast",
    "group", "ontology", "term", "module"
  ),
  definitions = c("module", "gene_symbol")
)

for (object_name in names(required_columns)) {
  object <- get(object_name)
  missing <- setdiff(required_columns[[object_name]], names(object))

  if (length(missing) > 0) {
    stop(
      "Colonne mancanti in ", object_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

integrated$gene_symbol <- toupper(clean_text(integrated$gene_symbol))
integrated$module <- clean_text(integrated$module)
overlap$gene_symbol <- toupper(clean_text(overlap$gene_symbol))
overlap$module <- clean_text(overlap$module)
definitions$gene_symbol <- toupper(clean_text(definitions$gene_symbol))
definitions$module <- clean_text(definitions$module)

numeric_columns <- c(
  "logFC_siNRP1_3", "FDR_siNRP1_3",
  "logFC_siNRP1_5", "FDR_siNRP1_5"
)

for (column in numeric_columns) {
  integrated[[column]] <- suppressWarnings(
    as.numeric(integrated[[column]])
  )
}

logical_columns <- c("significant_siNRP1_3", "significant_siNRP1_5")

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  value <- toupper(trimws(as.character(x)))
  result <- rep(NA, length(value))
  result[value %in% c("TRUE", "T", "1", "YES")] <- TRUE
  result[value %in% c("FALSE", "F", "0", "NO")] <- FALSE
  result
}

for (column in logical_columns) {
  integrated[[column]] <- to_logical(integrated[[column]])
}

# Rimuove solo righe completamente prive delle chiavi necessarie.
integrated <- integrated[
  !is.na(integrated$gene_symbol) & !is.na(integrated$module),
  ,
  drop = FALSE
]

overlap <- overlap[
  !is.na(overlap$gene_symbol) &
    !is.na(overlap$module) &
    !is.na(overlap$term),
  ,
  drop = FALSE
]

definitions <- unique(definitions[
  !is.na(definitions$gene_symbol) & !is.na(definitions$module),
  c("module", "gene_symbol"),
  drop = FALSE
])

# ------------------------------------------------------------
# 4. GENE-MODULE MEMBERSHIP AND BRIDGE GENES
# ------------------------------------------------------------

membership_split <- split(
  definitions$module,
  definitions$gene_symbol
)

membership <- data.frame(
  gene_symbol = names(membership_split),
  module_count = vapply(
    membership_split,
    function(x) length(unique(x)),
    integer(1)
  ),
  modules = vapply(
    membership_split,
    function(x) paste(sort(unique(x)), collapse = " | "),
    character(1)
  ),
  stringsAsFactors = FALSE
)

integrated_gene_level <- integrated[
  !duplicated(integrated$gene_symbol),
  ,
  drop = FALSE
]

bridge_genes <- merge(
  membership[membership$module_count > 1, , drop = FALSE],
  integrated_gene_level,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

bridge_genes <- bridge_genes[
  order(-bridge_genes$module_count, bridge_genes$gene_symbol),
  ,
  drop = FALSE
]

write.csv(
  membership,
  file.path(OUTPUT_TABLES, "gene_module_membership.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  bridge_genes,
  file.path(OUTPUT_TABLES, "cross_module_bridge_genes.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 5. INTEGRATE PATHWAY SIGNATURES WITH MODULE DE RESULTS
# ------------------------------------------------------------

# Chiave composta: un gene può appartenere legittimamente a più moduli.
integrated_key <- paste(
  integrated$gene_symbol,
  integrated$module,
  sep = "|||"
)

if (anyDuplicated(integrated_key)) {
  stop(
    "La tabella integrata contiene più righe per la stessa coppia ",
    "gene_symbol + module. Verificare l'output dello script 33."
  )
}

overlap_joined <- merge(
  overlap,
  integrated,
  by = c("gene_symbol", "module"),
  all.x = TRUE,
  sort = FALSE
)

overlap_joined$de_record_available <- !is.na(
  overlap_joined$logFC_siNRP1_3
) | !is.na(overlap_joined$logFC_siNRP1_5)

write.csv(
  overlap_joined,
  file.path(OUTPUT_TABLES, "pathway_module_gene_integration.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 6. PATHWAY x MODULE SUMMARY
# ------------------------------------------------------------

group_columns <- c(
  "signature_row", "contrast", "group",
  "ontology", "term", "module"
)

group_key <- do.call(
  interaction,
  c(
    unname(as.list(overlap_joined[, group_columns, drop = FALSE])),
    list(drop = TRUE, lex.order = TRUE, sep = "|||")
  )
)

split_rows <- split(
  seq_len(nrow(overlap_joined)),
  group_key
)

summarize_signature_module <- function(index) {
  x <- overlap_joined[index, , drop = FALSE]

  genes <- sort(unique(x$gene_symbol))
  sig3 <- !is.na(x$significant_siNRP1_3) &
    x$significant_siNRP1_3
  sig5 <- !is.na(x$significant_siNRP1_5) &
    x$significant_siNRP1_5

  both <- sig3 & sig5
  concordant <- both &
    x$direction_concordance == "CONCORDANT"
  discordant <- both &
    x$direction_concordance == "DISCORDANT"

  data.frame(
    signature_row = x$signature_row[1],
    contrast = x$contrast[1],
    group = x$group[1],
    ontology = x$ontology[1],
    term = x$term[1],
    module = x$module[1],
    signature_module_gene_count = length(genes),
    signature_module_genes = paste(genes, collapse = ";"),
    genes_detected_in_DE = length(unique(
      x$gene_symbol[x$de_record_available]
    )),
    significant_siNRP1_3 = length(unique(x$gene_symbol[sig3])),
    significant_siNRP1_5 = length(unique(x$gene_symbol[sig5])),
    significant_both = length(unique(x$gene_symbol[both])),
    concordant_significant_both = length(unique(
      x$gene_symbol[concordant]
    )),
    discordant_significant_both = length(unique(
      x$gene_symbol[discordant]
    )),
    up_siNRP1_3 = length(unique(
      x$gene_symbol[sig3 & x$direction_siNRP1_3 == "UP"]
    )),
    down_siNRP1_3 = length(unique(
      x$gene_symbol[sig3 & x$direction_siNRP1_3 == "DOWN"]
    )),
    up_siNRP1_5 = length(unique(
      x$gene_symbol[sig5 & x$direction_siNRP1_5 == "UP"]
    )),
    down_siNRP1_5 = length(unique(
      x$gene_symbol[sig5 & x$direction_siNRP1_5 == "DOWN"]
    )),
    stringsAsFactors = FALSE
  )
}

pathway_module_summary <- do.call(
  rbind,
  lapply(split_rows, summarize_signature_module)
)

rownames(pathway_module_summary) <- NULL

write.csv(
  pathway_module_summary,
  file.path(OUTPUT_TABLES, "pathway_module_support_summary.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 7. SHARED VS KNOCKDOWN-SPECIFIC PATHWAY-MODULE SIGNALS
# ------------------------------------------------------------

term_keys <- c("group", "ontology", "term", "module")

term_presence <- unique(pathway_module_summary[
  ,
  c(term_keys, "contrast"),
  drop = FALSE
])

contrast_3 <- "siNRP1_3_vs_siNT"
contrast_5 <- "siNRP1_5_vs_siNT"

present_3 <- term_presence[
  term_presence$contrast == contrast_3,
  term_keys,
  drop = FALSE
]
present_5 <- term_presence[
  term_presence$contrast == contrast_5,
  term_keys,
  drop = FALSE
]

present_3 <- unique(present_3)
present_5 <- unique(present_5)
present_3$present_siNRP1_3 <- TRUE
present_5$present_siNRP1_5 <- TRUE

pathway_comparison <- merge(
  present_3,
  present_5,
  by = term_keys,
  all = TRUE,
  sort = FALSE
)

pathway_comparison$present_siNRP1_3[
  is.na(pathway_comparison$present_siNRP1_3)
] <- FALSE
pathway_comparison$present_siNRP1_5[
  is.na(pathway_comparison$present_siNRP1_5)
] <- FALSE

pathway_comparison$comparison_class <- ifelse(
  pathway_comparison$present_siNRP1_3 &
    pathway_comparison$present_siNRP1_5,
  "SHARED",
  ifelse(
    pathway_comparison$present_siNRP1_3,
    "ONLY_siNRP1_3",
    "ONLY_siNRP1_5"
  )
)

write.csv(
  pathway_comparison,
  file.path(OUTPUT_TABLES, "shared_vs_specific_pathway_module_signals.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 8. MODULE-LEVEL SUMMARY OF PATHWAY SUPPORT
# ------------------------------------------------------------

module_pathway_counts <- aggregate(
  term ~ module + contrast + group + ontology,
  data = unique(pathway_module_summary[
    ,
    c("module", "contrast", "group", "ontology", "term"),
    drop = FALSE
  ]),
  FUN = function(x) length(unique(x))
)

names(module_pathway_counts)[
  names(module_pathway_counts) == "term"
] <- "distinct_pathway_terms"

write.csv(
  module_pathway_counts,
  file.path(OUTPUT_TABLES, "module_pathway_counts.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 9. ESSENTIAL AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "integrated_input_rows",
    "integrated_unique_genes",
    "integrated_gene_module_assignments",
    "module_definition_assignments",
    "signature_module_overlap_rows",
    "signature_module_summary_rows",
    "overlap_rows_with_DE_support",
    "unique_bridge_genes",
    "pathway_module_comparison_rows",
    "shared_pathway_module_signals",
    "siNRP1_3_specific_pathway_module_signals",
    "siNRP1_5_specific_pathway_module_signals"
  ),
  value = c(
    nrow(integrated),
    length(unique(integrated$gene_symbol)),
    nrow(integrated),
    nrow(definitions),
    nrow(overlap),
    nrow(pathway_module_summary),
    sum(overlap_joined$de_record_available, na.rm = TRUE),
    nrow(bridge_genes),
    nrow(pathway_comparison),
    sum(pathway_comparison$comparison_class == "SHARED"),
    sum(pathway_comparison$comparison_class == "ONLY_siNRP1_3"),
    sum(pathway_comparison$comparison_class == "ONLY_siNRP1_5")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(OUTPUT_AUDIT, "script34_audit.csv"),
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 10. TEXT REPORT
# ------------------------------------------------------------

report_path <- file.path(
  OUTPUT_REPORTS,
  "integrated_pathway_module_analysis_report.txt"
)

sink(report_path, split = TRUE)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("INTEGRATED PATHWAY-MODULE-DE ANALYSIS - SCRIPT 34\n")
cat("============================================================\n\n")

cat("PROJECT:", PROJECT_ROOT, "\n")
cat("OUTPUT:", OUTPUT_ROOT, "\n\n")

cat("AUDIT\n")
cat("------------------------------------------------------------\n")
print(audit, row.names = FALSE)

cat("\nPATHWAY-MODULE SIGNALS BY CONTRAST\n")
cat("------------------------------------------------------------\n")
print(
  aggregate(
    cbind(
      significant_siNRP1_3,
      significant_siNRP1_5,
      significant_both,
      concordant_significant_both,
      discordant_significant_both
    ) ~ contrast,
    data = pathway_module_summary,
    FUN = sum
  ),
  row.names = FALSE
)

cat("\nSHARED VS SPECIFIC SIGNALS\n")
cat("------------------------------------------------------------\n")
print(
  as.data.frame(table(pathway_comparison$comparison_class)),
  row.names = FALSE
)

cat("\nCROSS-MODULE BRIDGE GENES\n")
cat("------------------------------------------------------------\n")
if (nrow(bridge_genes) == 0) {
  cat("Nessun gene ponte rilevato.\n")
} else {
  print(
    bridge_genes[
      ,
      intersect(
        c(
          "gene_symbol", "module_count", "modules",
          "logFC_siNRP1_3", "FDR_siNRP1_3",
          "logFC_siNRP1_5", "FDR_siNRP1_5",
          "joint_class"
        ),
        names(bridge_genes)
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

cat("\nINTERPRETATION LIMITATIONS\n")
cat("------------------------------------------------------------\n")
cat("1. Le sovrapposizioni pathway-modulo sono descrittive.\n")
cat("2. Non sono stati eseguiti nuovi test di arricchimento.\n")
cat("3. La presenza in una signature non dimostra attivazione del pathway.\n")
cat("4. I geni possono appartenere a più moduli; i conteggi per modulo\n")
cat("   non devono essere sommati per ottenere geni unici.\n")
cat("5. La concordanza descrive il segno del logFC nei due contrasti,\n")
cat("   non dimostra un effetto causale comune.\n")

cat("\nOUTPUT FILES\n")
cat("------------------------------------------------------------\n")
for (directory in c(OUTPUT_TABLES, OUTPUT_AUDIT, OUTPUT_REPORTS)) {
  files <- list.files(directory, full.names = FALSE)
  cat("\n", directory, "\n", sep = "")
  cat(paste(files, collapse = "\n"), "\n")
}

sink()

cat("\nSCRIPT 34 COMPLETATO\n")
cat("Output:", OUTPUT_ROOT, "\n")
cat("Report:", report_path, "\n")
cat("Tabelle create:", length(list.files(OUTPUT_TABLES)), "\n")
cat("Audit creati:", length(list.files(OUTPUT_AUDIT)), "\n")
