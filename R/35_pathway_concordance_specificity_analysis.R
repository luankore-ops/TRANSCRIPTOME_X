
# ============================================================
# TRANSCRIPTOME-X
# Script 35
# PATHWAY CONCORDANCE AND SPECIFICITY ANALYSIS
#
# Confronto descrittivo delle signature funzionali:
# - pathway condivisi e specifici
# - supporto dei geni DE
# - concordanza / discordanza dei logFC
# - integrazione con i moduli biologici
#
# Nessun nuovo test di arricchimento.
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("PATHWAY CONCORDANCE AND SPECIFICITY - SCRIPT 35\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------

find_project_root <- function() {
  
  candidates <- c(
    getwd(),
    file.path(getwd(), "TRANSCRIPTOME_X"),
    "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"
  )
  
  for (candidate in unique(candidates)) {
    
    if (
      dir.exists(file.path(candidate, "R")) &&
      dir.exists(file.path(candidate, "results"))
    ) {
      return(normalizePath(
        candidate,
        winslash = "/",
        mustWork = TRUE
      ))
    }
  }
  
  stop("Directory principale del progetto non trovata.")
}

PROJECT <- find_project_root()

INPUT_SIGNATURES <- file.path(
  PROJECT,
  "results/functional_enrichment_interpretation/tables",
  "significant_term_gene_signatures.csv"
)

INPUT_DE_3 <- file.path(
  PROJECT,
  "results/normalization_sensitivity",
  "DE_TMM_siNRP1_3_vs_siNT.csv"
)

INPUT_DE_5 <- file.path(
  PROJECT,
  "results/normalization_sensitivity",
  "DE_TMM_siNRP1_5_vs_siNT.csv"
)

INPUT_MODULES <- file.path(
  PROJECT,
  "results/integrated_gene_module_analysis/tables",
  "module_gene_definitions.csv"
)

OUTPUT_ROOT <- file.path(
  PROJECT,
  "results/pathway_concordance_specificity"
)

OUTPUT_TABLES <- file.path(OUTPUT_ROOT, "tables")
OUTPUT_AUDIT <- file.path(OUTPUT_ROOT, "audit")
OUTPUT_REPORTS <- file.path(OUTPUT_ROOT, "reports")

dir.create(OUTPUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_AUDIT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_REPORTS, recursive = TRUE, showWarnings = FALSE)

INPUT_FILES <- c(
  signatures = INPUT_SIGNATURES,
  DE_siNRP1_3 = INPUT_DE_3,
  DE_siNRP1_5 = INPUT_DE_5,
  module_definitions = INPUT_MODULES
)

missing_files <- INPUT_FILES[!file.exists(INPUT_FILES)]

if (length(missing_files) > 0) {
  stop(
    "Input mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("PROJECT:", PROJECT, "\n")
cat("OUTPUT:", OUTPUT_ROOT, "\n\n")

# ------------------------------------------------------------
# 2. GENERAL FUNCTIONS
# ------------------------------------------------------------

read_csv_checked <- function(path) {
  
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

find_column <- function(data, candidates, description) {
  
  index <- match(
    tolower(candidates),
    tolower(names(data))
  )
  
  index <- index[!is.na(index)]
  
  if (length(index) == 0) {
    stop(
      "Colonna non trovata: ", description,
      "\nDisponibili: ",
      paste(names(data), collapse = ", ")
    )
  }
  
  names(data)[index[1]]
}

clean_text <- function(x) {
  
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

normalize_symbol <- function(x) {
  toupper(clean_text(x))
}

normalize_description <- function(x) {
  
  x <- tolower(clean_text(x))
  x <- gsub("[[:punct:]]+", " ", x)
  x <- gsub("\\s+", " ", x)
  
  trimws(x)
}

parse_gene_list <- function(x) {
  
  if (
    length(x) == 0 ||
    is.na(x) ||
    !nzchar(trimws(x))
  ) {
    return(character(0))
  }
  
  genes <- unlist(
    strsplit(as.character(x), "[,;/|]+")
  )
  
  genes <- normalize_symbol(genes)
  
  unique(genes[nzchar(genes)])
}

collapse_values <- function(x) {
  
  x <- unique(clean_text(x))
  x <- sort(x[nzchar(x)])
  
  paste(x, collapse = " | ")
}

write_output <- function(data, filename) {
  
  write.csv(
    data,
    file.path(OUTPUT_TABLES, filename),
    row.names = FALSE,
    na = ""
  )
  
  cat("Creato:", filename, "\n")
}

# ------------------------------------------------------------
# 3. LOAD INPUTS
# ------------------------------------------------------------

cat("CARICAMENTO INPUT\n")
cat("------------------------------------------------------------\n")

signatures_raw <- read_csv_checked(INPUT_SIGNATURES)
de3_raw <- read_csv_checked(INPUT_DE_3)
de5_raw <- read_csv_checked(INPUT_DE_5)
modules_raw <- read_csv_checked(INPUT_MODULES)

cat("Signature rows:", nrow(signatures_raw), "\n")
cat("DE siNRP1.3 rows:", nrow(de3_raw), "\n")
cat("DE siNRP1.5 rows:", nrow(de5_raw), "\n")
cat("Module definitions:", nrow(modules_raw), "\n\n")

# ------------------------------------------------------------
# 4. STANDARDIZE SIGNATURES
# ------------------------------------------------------------

gene_col <- find_column(
  signatures_raw,
  c("genes", "geneID", "gene_ids", "gene_symbols"),
  "signature genes"
)

term_col <- find_column(
  signatures_raw,
  c("Description", "description", "term"),
  "pathway description"
)

contrast_col <- find_column(
  signatures_raw,
  c("contrast"),
  "contrast"
)

group_col <- find_column(
  signatures_raw,
  c("group"),
  "group"
)

ontology_col <- find_column(
  signatures_raw,
  c("ontology"),
  "ontology"
)

signatures <- data.frame(
  genes_raw = clean_text(signatures_raw[[gene_col]]),
  term = clean_text(signatures_raw[[term_col]]),
  contrast_raw = clean_text(signatures_raw[[contrast_col]]),
  group = clean_text(signatures_raw[[group_col]]),
  ontology = clean_text(signatures_raw[[ontology_col]]),
  stringsAsFactors = FALSE
)

signatures$term_normalized <- normalize_description(
  signatures$term
)

classify_contrast <- function(x) {
  
  x <- toupper(gsub("[^A-Z0-9]", "", x))
  
  if (grepl("NRP13", x, fixed = TRUE)) {
    return("siNRP1_3")
  }
  
  if (grepl("NRP15", x, fixed = TRUE)) {
    return("siNRP1_5")
  }
  
  NA_character_
}

signatures$contrast <- vapply(
  signatures$contrast_raw,
  classify_contrast,
  character(1)
)

if (anyNA(signatures$contrast)) {
  
  stop(
    "Contrasti non riconosciuti: ",
    paste(
      unique(signatures$contrast_raw[is.na(signatures$contrast)]),
      collapse = ", "
    )
  )
}

signatures <- signatures[
  nzchar(signatures$term) &
    nzchar(signatures$group) &
    nzchar(signatures$ontology) &
    nzchar(signatures$term_normalized),
  ,
  drop = FALSE
]

# Identità pathway:
# ontologia + gruppo + descrizione normalizzata.

signatures$pathway_id <- paste(
  toupper(signatures$ontology),
  toupper(signatures$group),
  signatures$term_normalized,
  sep = "|||"
)

# ------------------------------------------------------------
# 5. EXPAND SIGNATURES TO GENE LEVEL
# ------------------------------------------------------------

signature_list <- vector("list", nrow(signatures))

for (i in seq_len(nrow(signatures))) {
  
  genes <- parse_gene_list(signatures$genes_raw[i])
  
  if (length(genes) == 0) {
    next
  }
  
  signature_list[[i]] <- data.frame(
    pathway_id = signatures$pathway_id[i],
    term = signatures$term[i],
    ontology = signatures$ontology[i],
    group = signatures$group[i],
    contrast = signatures$contrast[i],
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )
}

signature_list <- Filter(
  Negate(is.null),
  signature_list
)

if (length(signature_list) == 0) {
  stop("Nessun gene trovato nelle signature.")
}

signature_genes <- unique(
  do.call(rbind, signature_list)
)

pathway_metadata <- unique(
  signatures[
    ,
    c(
      "pathway_id",
      "term",
      "term_normalized",
      "ontology",
      "group"
    ),
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 6. PREPARE DIFFERENTIAL EXPRESSION
# ------------------------------------------------------------

prepare_de <- function(data, contrast_name) {
  
  symbol_col <- find_column(
    data,
    c(
      "Ensembl.114.Gene.Name",
      "gene_symbol",
      "symbol",
      "gene_name"
    ),
    paste(contrast_name, "gene symbol")
  )
  
  logfc_col <- find_column(
    data,
    c("logFC", "log2FoldChange"),
    paste(contrast_name, "logFC")
  )
  
  fdr_col <- find_column(
    data,
    c("FDR", "padj", "adjusted_p_value"),
    paste(contrast_name, "FDR")
  )
  
  result <- data.frame(
    gene_symbol = normalize_symbol(data[[symbol_col]]),
    logFC = suppressWarnings(
      as.numeric(data[[logfc_col]])
    ),
    FDR = suppressWarnings(
      as.numeric(data[[fdr_col]])
    ),
    stringsAsFactors = FALSE
  )
  
  result <- result[
    nzchar(result$gene_symbol) &
      is.finite(result$logFC),
    ,
    drop = FALSE
  ]
  
  result$FDR_order <- result$FDR
  result$FDR_order[!is.finite(result$FDR_order)] <- Inf
  
  result <- result[
    order(result$gene_symbol, result$FDR_order),
    ,
    drop = FALSE
  ]
  
  result <- result[
    !duplicated(result$gene_symbol),
    c("gene_symbol", "logFC", "FDR"),
    drop = FALSE
  ]
  
  names(result)[2:3] <- c(
    paste0("logFC_", contrast_name),
    paste0("FDR_", contrast_name)
  )
  
  fdr_name <- paste0("FDR_", contrast_name)
  
  result[[paste0("significant_", contrast_name)]] <-
    !is.na(result[[fdr_name]]) &
    result[[fdr_name]] < 0.05
  
  result
}

de3 <- prepare_de(de3_raw, "siNRP1_3")
de5 <- prepare_de(de5_raw, "siNRP1_5")

cat("Geni DE unici siNRP1.3:", nrow(de3), "\n")
cat("Geni DE unici siNRP1.5:", nrow(de5), "\n\n")

# ------------------------------------------------------------
# 7. PREPARE MODULE DEFINITIONS
# ------------------------------------------------------------

module_gene_col <- find_column(
  modules_raw,
  c(
    "gene_symbol",
    "gene",
    "symbol",
    "Ensembl.114.Gene.Name"
  ),
  "module gene symbol"
)

module_name_col <- find_column(
  modules_raw,
  c("module", "module_name"),
  "module name"
)

modules <- unique(
  data.frame(
    gene_symbol = normalize_symbol(
      modules_raw[[module_gene_col]]
    ),
    module = clean_text(
      modules_raw[[module_name_col]]
    ),
    stringsAsFactors = FALSE
  )
)

modules <- modules[
  nzchar(modules$gene_symbol) &
    nzchar(modules$module),
  ,
  drop = FALSE
]

module_membership <- aggregate(
  module ~ gene_symbol,
  data = modules,
  FUN = collapse_values
)

names(module_membership)[2] <- "modules"

# ------------------------------------------------------------
# 8. BUILD PATHWAY-GENE UNIVERSE
# ------------------------------------------------------------

pathway_gene_universe <- unique(
  signature_genes[
    ,
    c("pathway_id", "gene_symbol"),
    drop = FALSE
  ]
)

presence3 <- unique(
  signature_genes[
    signature_genes$contrast == "siNRP1_3",
    c("pathway_id", "gene_symbol"),
    drop = FALSE
  ]
)

presence5 <- unique(
  signature_genes[
    signature_genes$contrast == "siNRP1_5",
    c("pathway_id", "gene_symbol"),
    drop = FALSE
  ]
)

presence3$in_signature_siNRP1_3 <- TRUE
presence5$in_signature_siNRP1_5 <- TRUE

pathway_gene_evidence <- merge(
  pathway_gene_universe,
  pathway_metadata,
  by = "pathway_id",
  all.x = TRUE
)

pathway_gene_evidence <- merge(
  pathway_gene_evidence,
  de3,
  by = "gene_symbol",
  all.x = TRUE
)

pathway_gene_evidence <- merge(
  pathway_gene_evidence,
  de5,
  by = "gene_symbol",
  all.x = TRUE
)

pathway_gene_evidence <- merge(
  pathway_gene_evidence,
  presence3,
  by = c("pathway_id", "gene_symbol"),
  all.x = TRUE
)

pathway_gene_evidence <- merge(
  pathway_gene_evidence,
  presence5,
  by = c("pathway_id", "gene_symbol"),
  all.x = TRUE
)

for (column in c(
  "in_signature_siNRP1_3",
  "in_signature_siNRP1_5",
  "significant_siNRP1_3",
  "significant_siNRP1_5"
)) {
  
  pathway_gene_evidence[[column]][
    is.na(pathway_gene_evidence[[column]])
  ] <- FALSE
}

# ------------------------------------------------------------
# 9. GENE-LEVEL CONCORDANCE
# ------------------------------------------------------------

pathway_gene_evidence$significant_both <-
  pathway_gene_evidence$significant_siNRP1_3 &
  pathway_gene_evidence$significant_siNRP1_5

pathway_gene_evidence$direction_siNRP1_3 <- ifelse(
  is.na(pathway_gene_evidence$logFC_siNRP1_3),
  NA_character_,
  ifelse(
    pathway_gene_evidence$logFC_siNRP1_3 > 0,
    "UP",
    ifelse(
      pathway_gene_evidence$logFC_siNRP1_3 < 0,
      "DOWN",
      "ZERO"
    )
  )
)

pathway_gene_evidence$direction_siNRP1_5 <- ifelse(
  is.na(pathway_gene_evidence$logFC_siNRP1_5),
  NA_character_,
  ifelse(
    pathway_gene_evidence$logFC_siNRP1_5 > 0,
    "UP",
    ifelse(
      pathway_gene_evidence$logFC_siNRP1_5 < 0,
      "DOWN",
      "ZERO"
    )
  )
)

pathway_gene_evidence$joint_class <- "NOT_SIGNIFICANT_BOTH"

both <- pathway_gene_evidence$significant_both

pathway_gene_evidence$joint_class[
  both &
    pathway_gene_evidence$direction_siNRP1_3 ==
    pathway_gene_evidence$direction_siNRP1_5
] <- "SIGNIFICANT_BOTH_CONCORDANT"

pathway_gene_evidence$joint_class[
  both &
    pathway_gene_evidence$direction_siNRP1_3 !=
    pathway_gene_evidence$direction_siNRP1_5
] <- "SIGNIFICANT_BOTH_DISCORDANT"

pathway_gene_evidence$joint_class[
  pathway_gene_evidence$significant_siNRP1_3 &
    !pathway_gene_evidence$significant_siNRP1_5
] <- "SIGNIFICANT_ONLY_siNRP1_3"

pathway_gene_evidence$joint_class[
  !pathway_gene_evidence$significant_siNRP1_3 &
    pathway_gene_evidence$significant_siNRP1_5
] <- "SIGNIFICANT_ONLY_siNRP1_5"

# ------------------------------------------------------------
# 10. MODULE MEMBERSHIP
# ------------------------------------------------------------

pathway_gene_evidence <- merge(
  pathway_gene_evidence,
  module_membership,
  by = "gene_symbol",
  all.x = TRUE
)

pathway_gene_evidence$modules[
  is.na(pathway_gene_evidence$modules)
] <- "NO_DEFINED_MODULE"

# ------------------------------------------------------------
# 11. PATHWAY-LEVEL SUMMARY
# ------------------------------------------------------------

summarize_pathway <- function(pathway_id) {
  
  d <- pathway_gene_evidence[
    pathway_gene_evidence$pathway_id == pathway_id,
    ,
    drop = FALSE
  ]
  
  meta <- pathway_metadata[
    pathway_metadata$pathway_id == pathway_id,
    ,
    drop = FALSE
  ]
  
  sig3 <- unique(
    d$gene_symbol[d$in_signature_siNRP1_3]
  )
  
  sig5 <- unique(
    d$gene_symbol[d$in_signature_siNRP1_5]
  )
  
  de3_genes <- unique(
    d$gene_symbol[d$significant_siNRP1_3]
  )
  
  de5_genes <- unique(
    d$gene_symbol[d$significant_siNRP1_5]
  )
  
  both_genes <- unique(
    d$gene_symbol[d$significant_both]
  )
  
  concordant <- unique(
    d$gene_symbol[
      d$joint_class == "SIGNIFICANT_BOTH_CONCORDANT"
    ]
  )
  
  discordant <- unique(
    d$gene_symbol[
      d$joint_class == "SIGNIFICANT_BOTH_DISCORDANT"
    ]
  )
  
  data.frame(
    pathway_id = pathway_id,
    term = meta$term[1],
    ontology = meta$ontology[1],
    group = meta$group[1],
    
    present_siNRP1_3 = length(sig3) > 0,
    present_siNRP1_5 = length(sig5) > 0,
    
    signature_genes_siNRP1_3 = length(sig3),
    signature_genes_siNRP1_5 = length(sig5),
    
    unique_signature_genes = length(
      unique(c(sig3, sig5))
    ),
    
    significant_genes_siNRP1_3 = length(de3_genes),
    significant_genes_siNRP1_5 = length(de5_genes),
    
    significant_both = length(both_genes),
    concordant_significant_both = length(concordant),
    discordant_significant_both = length(discordant),
    
    concordance_fraction = ifelse(
      length(both_genes) > 0,
      length(concordant) / length(both_genes),
      NA_real_
    ),
    
    stringsAsFactors = FALSE
  )
}

pathway_summary <- do.call(
  rbind,
  lapply(
    unique(pathway_metadata$pathway_id),
    summarize_pathway
  )
)

pathway_summary$pathway_status <- ifelse(
  pathway_summary$present_siNRP1_3 &
    pathway_summary$present_siNRP1_5,
  "SHARED",
  ifelse(
    pathway_summary$present_siNRP1_3,
    "ONLY_siNRP1_3",
    "ONLY_siNRP1_5"
  )
)

pathway_summary <- pathway_summary[
  order(
    pathway_summary$pathway_status,
    pathway_summary$ontology,
    pathway_summary$group,
    pathway_summary$term
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 12. SHARED VS SPECIFIC PATHWAYS
# ------------------------------------------------------------

shared_specificity <- pathway_summary[
  ,
  c(
    "pathway_id",
    "term",
    "ontology",
    "group",
    "pathway_status",
    "signature_genes_siNRP1_3",
    "signature_genes_siNRP1_5",
    "unique_signature_genes",
    "significant_genes_siNRP1_3",
    "significant_genes_siNRP1_5",
    "significant_both",
    "concordant_significant_both",
    "discordant_significant_both",
    "concordance_fraction"
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# 13. PATHWAY-MODULE CROSSWALK
# ------------------------------------------------------------

module_records <- pathway_gene_evidence[
  pathway_gene_evidence$modules != "NO_DEFINED_MODULE",
  ,
  drop = FALSE
]

module_crosswalk <- data.frame()

if (nrow(module_records) > 0) {
  
  module_list <- strsplit(
    module_records$modules,
    " \\| ",
    perl = TRUE
  )
  
  expanded_index <- rep(
    seq_len(nrow(module_records)),
    lengths(module_list)
  )
  
  expanded_modules <- unlist(
    module_list,
    use.names = FALSE
  )
  
  expanded <- module_records[
    expanded_index,
    ,
    drop = FALSE
  ]
  
  expanded$module <- expanded_modules
  
  expanded <- unique(
    expanded[
      ,
      c(
        "pathway_id",
        "term",
        "ontology",
        "group",
        "module",
        "gene_symbol",
        "significant_siNRP1_3",
        "significant_siNRP1_5",
        "significant_both",
        "joint_class"
      ),
      drop = FALSE
    ]
  )
  
  module_crosswalk <- aggregate(
    cbind(
      gene_count = rep(1, nrow(expanded)),
      significant_siNRP1_3 = as.numeric(
        expanded$significant_siNRP1_3
      ),
      significant_siNRP1_5 = as.numeric(
        expanded$significant_siNRP1_5
      ),
      significant_both = as.numeric(
        expanded$significant_both
      ),
      concordant = as.numeric(
        expanded$joint_class ==
          "SIGNIFICANT_BOTH_CONCORDANT"
      ),
      discordant = as.numeric(
        expanded$joint_class ==
          "SIGNIFICANT_BOTH_DISCORDANT"
      )
    ),
    by = list(
      pathway_id = expanded$pathway_id,
      term = expanded$term,
      ontology = expanded$ontology,
      group = expanded$group,
      module = expanded$module
    ),
    FUN = sum
  )
}

# ------------------------------------------------------------
# 14. PATHWAY GENE LISTS
# ------------------------------------------------------------

pathway_gene_lists <- aggregate(
  gene_symbol ~ pathway_id,
  data = pathway_gene_evidence,
  FUN = collapse_values
)

names(pathway_gene_lists)[2] <- "pathway_genes"

pathway_gene_lists <- merge(
  pathway_metadata,
  pathway_gene_lists,
  by = "pathway_id",
  all.x = TRUE
)

# ------------------------------------------------------------
# 15. AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "signature_input_rows",
    "valid_signature_rows",
    "signature_gene_records",
    "unique_signature_genes",
    "unique_pathways",
    "pathways_siNRP1_3",
    "pathways_siNRP1_5",
    "shared_pathways",
    "siNRP1_3_specific_pathways",
    "siNRP1_5_specific_pathways",
    "pathway_gene_evidence_rows",
    "unique_genes_in_pathway_evidence",
    "significant_pathway_gene_records_both",
    "concordant_pathway_gene_records",
    "discordant_pathway_gene_records",
    "module_crosswalk_rows",
    "module_definition_rows"
  ),
  value = c(
    nrow(signatures_raw),
    nrow(signatures),
    nrow(signature_genes),
    length(unique(signature_genes$gene_symbol)),
    nrow(pathway_summary),
    sum(pathway_summary$present_siNRP1_3),
    sum(pathway_summary$present_siNRP1_5),
    sum(pathway_summary$pathway_status == "SHARED"),
    sum(pathway_summary$pathway_status == "ONLY_siNRP1_3"),
    sum(pathway_summary$pathway_status == "ONLY_siNRP1_5"),
    nrow(pathway_gene_evidence),
    length(unique(pathway_gene_evidence$gene_symbol)),
    sum(pathway_gene_evidence$significant_both),
    sum(
      pathway_gene_evidence$joint_class ==
        "SIGNIFICANT_BOTH_CONCORDANT"
    ),
    sum(
      pathway_gene_evidence$joint_class ==
        "SIGNIFICANT_BOTH_DISCORDANT"
    ),
    nrow(module_crosswalk),
    nrow(modules)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 16. EXPORT TABLES
# ------------------------------------------------------------

cat("\n")
cat("ESPORTAZIONE RISULTATI\n")
cat("------------------------------------------------------------\n")

write_output(
  pathway_gene_evidence,
  "pathway_gene_evidence.csv"
)

write_output(
  pathway_summary,
  "pathway_concordance_summary.csv"
)

write_output(
  shared_specificity,
  "pathway_shared_specificity.csv"
)

write_output(
  pathway_gene_lists,
  "pathway_gene_lists.csv"
)

write_output(
  module_crosswalk,
  "pathway_module_crosswalk.csv"
)

write.csv(
  audit,
  file.path(OUTPUT_AUDIT, "script35_audit.csv"),
  row.names = FALSE,
  na = ""
)

cat("Creato: script35_audit.csv\n")

# ------------------------------------------------------------
# 17. REPORT
# ------------------------------------------------------------

report_path <- file.path(
  OUTPUT_REPORTS,
  "pathway_concordance_specificity_report.txt"
)

status_table <- as.data.frame(
  table(pathway_summary$pathway_status)
)

names(status_table) <- c(
  "pathway_status",
  "count"
)

group_table <- aggregate(
  pathway_id ~ ontology + group + pathway_status,
  data = pathway_summary,
  FUN = function(x) length(unique(x))
)

names(group_table)[4] <- "pathway_count"

significant_pathways <- pathway_summary[
  pathway_summary$significant_both > 0,
  ,
  drop = FALSE
]

significant_pathways <- significant_pathways[
  order(
    -significant_pathways$concordant_significant_both,
    -significant_pathways$significant_both
  ),
  ,
  drop = FALSE
]

top_pathways <- head(
  significant_pathways[
    ,
    c(
      "term",
      "ontology",
      "group",
      "pathway_status",
      "significant_both",
      "concordant_significant_both",
      "discordant_significant_both",
      "concordance_fraction"
    ),
    drop = FALSE
  ],
  20
)

gene_metrics <- data.frame(
  metric = c(
    "Pathway-gene records",
    "Unique genes",
    "Significant siNRP1.3 records",
    "Significant siNRP1.5 records",
    "Significant in both",
    "Concordant records",
    "Discordant records"
  ),
  value = c(
    nrow(pathway_gene_evidence),
    length(unique(pathway_gene_evidence$gene_symbol)),
    sum(pathway_gene_evidence$significant_siNRP1_3),
    sum(pathway_gene_evidence$significant_siNRP1_5),
    sum(pathway_gene_evidence$significant_both),
    sum(
      pathway_gene_evidence$joint_class ==
        "SIGNIFICANT_BOTH_CONCORDANT"
    ),
    sum(
      pathway_gene_evidence$joint_class ==
        "SIGNIFICANT_BOTH_DISCORDANT"
    )
  ),
  stringsAsFactors = FALSE
)

report_lines <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "PATHWAY CONCORDANCE AND SPECIFICITY - SCRIPT 35",
  "============================================================",
  "",
  paste("PROJECT:", PROJECT),
  paste("OUTPUT:", OUTPUT_ROOT),
  "",
  "PATHWAY STATUS",
  "------------------------------------------------------------",
  capture.output(print(status_table, row.names = FALSE)),
  "",
  "PATHWAY COUNTS BY ONTOLOGY AND GROUP",
  "------------------------------------------------------------",
  capture.output(print(group_table, row.names = FALSE)),
  "",
  "GENE-LEVEL SUPPORT",
  "------------------------------------------------------------",
  capture.output(print(gene_metrics, row.names = FALSE)),
  "",
  "TOP PATHWAYS WITH GENES SIGNIFICANT IN BOTH CONTRASTS",
  "------------------------------------------------------------",
  capture.output(print(top_pathways, row.names = FALSE)),
  "",
  "INTERPRETATION LIMITATIONS",
  "------------------------------------------------------------",
  "1. Il confronto usa le signature funzionali gia' disponibili.",
  "2. Non sono stati eseguiti nuovi test di arricchimento.",
  "3. SHARED indica la presenza dello stesso termine, ontologia",
  "   e gruppo nelle signature dei due confronti.",
  "4. ONLY indica assenza dalla signature dell'altro confronto,",
  "   non assenza biologica del pathway.",
  "5. I geni possono comparire in piu' pathway e moduli.",
  "6. I conteggi pathway-gene non sono conteggi di geni unici globali.",
  "7. Concordanza indica lo stesso segno del logFC nei due confronti.",
  "8. La significativita' e' definita come FDR < 0.05.",
  "9. L'appartenenza a una signature non dimostra attivazione.",
  "10. I risultati sono descrittivi e non dimostrano causalita'.",
  "",
  "OUTPUT FILES",
  "------------------------------------------------------------",
  "tables/pathway_gene_evidence.csv",
  "tables/pathway_concordance_summary.csv",
  "tables/pathway_shared_specificity.csv",
  "tables/pathway_gene_lists.csv",
  "tables/pathway_module_crosswalk.csv",
  "audit/script35_audit.csv",
  "reports/pathway_concordance_specificity_report.txt",
  "",
  "SCRIPT 35 COMPLETATO"
)

writeLines(report_lines, report_path)

# ------------------------------------------------------------
# 18. CONSOLE SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SCRIPT 35 - RISULTATI FINALI\n")
cat("============================================================\n\n")

cat("PATHWAY STATUS\n")
cat("------------------------------------------------------------\n")

print(status_table, row.names = FALSE)

cat("\nGENE-LEVEL SUPPORT\n")
cat("------------------------------------------------------------\n")

cat(
  "Pathway-gene records:",
  nrow(pathway_gene_evidence),
  "\n"
)

cat(
  "Unique genes:",
  length(unique(pathway_gene_evidence$gene_symbol)),
  "\n"
)

cat(
  "Significant in both:",
  sum(pathway_gene_evidence$significant_both),
  "\n"
)

cat(
  "Concordant:",
  sum(
    pathway_gene_evidence$joint_class ==
      "SIGNIFICANT_BOTH_CONCORDANT"
  ),
  "\n"
)

cat(
  "Discordant:",
  sum(
    pathway_gene_evidence$joint_class ==
      "SIGNIFICANT_BOTH_DISCORDANT"
  ),
  "\n"
)

cat("\nAUDIT\n")
cat("------------------------------------------------------------\n")

print(audit, row.names = FALSE)

cat("\n")
cat("============================================================\n")
cat("SCRIPT 35 COMPLETATO\n")
cat("Output:", OUTPUT_ROOT, "\n")
cat("Report:", report_path, "\n")
cat("Tabelle create: 5\n")
cat("Audit creati: 1\n")
cat("============================================================\n")
