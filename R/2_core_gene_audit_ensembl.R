
# ============================================================
# TRANSCRIPTOME-X - CORE GENE AUDIT
# Script: 22_core_gene_audit_ensembl.R
#
# Versione autonoma con validazione separata:
# - record Ensembl strict
# - record unresolved
#
# Non modifica gli output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE AUDIT V3\n")
cat("STRICT AND UNRESOLVED VALIDATION\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

input_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_analysis_v2"
)

output_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_audit_v3"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

integrated_file <- file.path(
  input_dir,
  "core_genes_with_differential_expression.csv"
)

membership_file <- file.path(
  input_dir,
  "gene_pathway_membership_ensembl.csv"
)

mapping_file <- file.path(
  input_dir,
  "gene_symbol_mapping_audit.csv"
)

strict_mapping_file <- file.path(
  input_dir,
  "gene_symbol_mapping_strict.csv"
)

required_files <- c(
  integrated_file,
  membership_file,
  mapping_file,
  strict_mapping_file
)

if (!all(file.exists(required_files))) {
  stop(
    "File di input mancanti:\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 2. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE\n")
cat("------------------------------------------------------------\n")

integrated <- read.csv(
  integrated_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

membership <- read.csv(
  membership_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping_audit <- read.csv(
  mapping_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

strict_mapping <- read.csv(
  strict_mapping_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Integrated:", nrow(integrated), "\n")
cat("Membership:", nrow(membership), "\n")
cat("Mapping audit:", nrow(mapping_audit), "\n")
cat("Strict mapping:", nrow(strict_mapping), "\n")

# ------------------------------------------------------------
# 3. VALIDAZIONE COLONNE
# ------------------------------------------------------------

required_integrated <- c(
  "comparison",
  "direction",
  "enrichment_symbol",
  "ensembl_id",
  "de_gene_name",
  "logFC",
  "FDR",
  "de_significant",
  "de_direction",
  "direction_consistent",
  "mapping_status",
  "de_mapping_status"
)

missing_columns <- setdiff(
  required_integrated,
  names(integrated)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti in integrated: ",
    paste(missing_columns, collapse = ", ")
  )
}

required_mapping <- c(
  "SYMBOL",
  "ENSEMBL"
)

missing_mapping <- setdiff(
  required_mapping,
  names(strict_mapping)
)

if (length(missing_mapping) > 0) {
  stop(
    "Colonne mancanti in strict mapping: ",
    paste(missing_mapping, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. NORMALIZZAZIONE ENSEMBL
# ------------------------------------------------------------

normalize_ensembl <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[is.na(x) | !nzchar(x)] <- NA_character_
  
  sub("\\.[0-9]+$", "", x)
}

integrated$ensembl_id <- normalize_ensembl(
  integrated$ensembl_id
)

strict_mapping$ENSEMBL <- normalize_ensembl(
  strict_mapping$ENSEMBL
)

# ------------------------------------------------------------
# 5. CLASSIFICAZIONE RECORD
# ------------------------------------------------------------

cat("\n2. CLASSIFICAZIONE RECORD\n")
cat("------------------------------------------------------------\n")

has_ensembl <- (
  !is.na(integrated$ensembl_id) &
    nzchar(integrated$ensembl_id)
)

strict_status <- (
  integrated$mapping_status == "STRICT_ONE_TO_ONE"
)

matched_status <- (
  integrated$de_mapping_status == "MATCHED_BY_ENSEMBL"
)

integrated$record_class <- ifelse(
  has_ensembl & strict_status & matched_status,
  "STRICT_INTEGRATED",
  "UNRESOLVED"
)

strict_rows <- integrated[
  integrated$record_class == "STRICT_INTEGRATED",
  ,
  drop = FALSE
]

unresolved_rows <- integrated[
  integrated$record_class == "UNRESOLVED",
  ,
  drop = FALSE
]

cat(
  "Record strict:",
  nrow(strict_rows),
  "\n"
)

cat(
  "Record unresolved:",
  nrow(unresolved_rows),
  "\n"
)

cat(
  "Totale:",
  nrow(integrated),
  "\n"
)

write.csv(
  integrated,
  file.path(output_dir, "all_core_gene_records_classified.csv"),
  row.names = FALSE
)

write.csv(
  strict_rows,
  file.path(output_dir, "strict_integrated_records.csv"),
  row.names = FALSE
)

write.csv(
  unresolved_rows,
  file.path(output_dir, "unresolved_records.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. CONTROLLO DUPLICATI
# ------------------------------------------------------------

cat("\n3. CONTROLLO DUPLICATI ENSEMBL\n")
cat("------------------------------------------------------------\n")

strict_key <- paste(
  strict_rows$comparison,
  strict_rows$direction,
  strict_rows$ensembl_id,
  sep = "|||"
)

duplicate_mask <- (
  duplicated(strict_key) |
    duplicated(strict_key, fromLast = TRUE)
)

duplicate_rows <- strict_rows[
  duplicate_mask,
  ,
  drop = FALSE
]

write.csv(
  duplicate_rows,
  file.path(output_dir, "duplicate_ensembl_core_gene_rows.csv"),
  row.names = FALSE
)

cat(
  "Duplicati:",
  nrow(duplicate_rows),
  "\n"
)

# ------------------------------------------------------------
# 7. AUDIT DE E DIREZIONE
# ------------------------------------------------------------

cat("\n4. AUDIT DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")

strict_rows$expected_direction <- ifelse(
  strict_rows$logFC > 0,
  "UP",
  ifelse(strict_rows$logFC < 0, "DOWN", "ZERO")
)

strict_rows$logFC_direction_match <- (
  strict_rows$expected_direction ==
    strict_rows$direction
)

strict_rows$DE_direction_match <- (
  strict_rows$de_direction ==
    strict_rows$direction
)

strict_rows$DE_significant_check <- (
  !is.na(strict_rows$de_significant) &
    strict_rows$de_significant == TRUE &
    !is.na(strict_rows$FDR) &
    strict_rows$FDR < 0.05
)

strict_rows$direction_consistency_check <- (
  strict_rows$logFC_direction_match &
    strict_rows$DE_direction_match &
    strict_rows$DE_significant_check
)

write.csv(
  strict_rows,
  file.path(output_dir, "strict_direction_consistency_audit.csv"),
  row.names = FALSE
)

unresolved_audit <- unresolved_rows

unresolved_audit$unresolved_reason <- ifelse(
  !has_ensembl[
    integrated$record_class == "UNRESOLVED"
  ],
  "MISSING_ENSEMBL",
  ifelse(
    !strict_status[
      integrated$record_class == "UNRESOLVED"
    ],
    "NON_STRICT_MAPPING",
    "DE_INTEGRATION_UNRESOLVED"
  )
)

write.csv(
  unresolved_audit,
  file.path(output_dir, "unresolved_mapping_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. VERIFICA UNICITA MAPPING STRICT
# ------------------------------------------------------------

cat("\n5. AUDIT MAPPING STRICT\n")
cat("------------------------------------------------------------\n")

strict_mapping <- strict_mapping[
  !is.na(strict_mapping$SYMBOL) &
    nzchar(strict_mapping$SYMBOL) &
    !is.na(strict_mapping$ENSEMBL) &
    nzchar(strict_mapping$ENSEMBL),
  ,
  drop = FALSE
]

symbol_counts <- table(strict_mapping$SYMBOL)
ensembl_counts <- table(strict_mapping$ENSEMBL)

symbol_unique <- all(symbol_counts == 1)
ensembl_unique <- all(ensembl_counts == 1)

write.csv(
  data.frame(
    SYMBOL = names(symbol_counts),
    ENSEMBL_count = as.integer(symbol_counts),
    unique_mapping = as.integer(symbol_counts) == 1
  ),
  file.path(output_dir, "strict_symbol_uniqueness.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    ENSEMBL = names(ensembl_counts),
    SYMBOL_count = as.integer(ensembl_counts),
    unique_mapping = as.integer(ensembl_counts) == 1
  ),
  file.path(output_dir, "strict_ensembl_uniqueness.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. AUDIT RECUPERO TRE GENI
# ------------------------------------------------------------

cat("\n6. AUDIT TRE GENI RECUPERATI\n")
cat("------------------------------------------------------------\n")

targets <- data.frame(
  comparison = c(
    "siNRP1_3_vs_siNT",
    "siNRP1_3_vs_siNT",
    "siNRP1_5_vs_siNT"
  ),
  enrichment_symbol = c(
    "SACK1G",
    "CDC50A",
    "CHLSN"
  ),
  expected_de_gene_name = c(
    "FAM83G",
    "TMEM30A",
    "C7orf50"
  ),
  expected_ensembl = c(
    "ENSG00000188522",
    "ENSG00000112697",
    "ENSG00000146540"
  ),
  stringsAsFactors = FALSE
)

recovery_audit <- merge(
  targets,
  integrated,
  by = c("comparison", "enrichment_symbol"),
  all.x = TRUE,
  sort = FALSE
)

recovery_audit$recovery_validated <- (
  recovery_audit$de_gene_name ==
    recovery_audit$expected_de_gene_name &
    recovery_audit$ensembl_id ==
    recovery_audit$expected_ensembl &
    recovery_audit$de_significant == TRUE &
    recovery_audit$de_mapping_status ==
    "MATCHED_BY_ENSEMBL"
)

write.csv(
  recovery_audit,
  file.path(output_dir, "three_gene_recovery_audit.csv"),
  row.names = FALSE
)

print(
  recovery_audit[
    ,
    c(
      "comparison",
      "enrichment_symbol",
      "de_gene_name",
      "ensembl_id",
      "logFC",
      "FDR",
      "recovery_validated"
    )
  ]
)

# ------------------------------------------------------------
# 10. RIEPILOGO PER CONFRONTO E DIREZIONE
# ------------------------------------------------------------

cat("\n7. RIEPILOGO PER CONFRONTO E DIREZIONE\n")
cat("------------------------------------------------------------\n")

summary_groups <- split(
  seq_len(nrow(integrated)),
  paste(
    integrated$comparison,
    integrated$direction,
    sep = "|||"
  )
)

audit_summary <- do.call(
  rbind,
  lapply(summary_groups, function(idx) {
    
    x <- integrated[idx, , drop = FALSE]
    
    strict_x <- x[
      x$record_class == "STRICT_INTEGRATED",
      ,
      drop = FALSE
    ]
    
    unresolved_x <- x[
      x$record_class == "UNRESOLVED",
      ,
      drop = FALSE
    ]
    
    data.frame(
      comparison = x$comparison[1],
      direction = x$direction[1],
      total_records = nrow(x),
      strict_records = nrow(strict_x),
      unresolved_records = nrow(unresolved_x),
      significant_DE = sum(
        strict_x$de_significant == TRUE,
        na.rm = TRUE
      ),
      direction_consistent = sum(
        strict_x$direction_consistency_check,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
)

rownames(audit_summary) <- NULL

print(audit_summary)

write.csv(
  audit_summary,
  file.path(output_dir, "core_gene_audit_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. VALIDAZIONE GLOBALE
# ------------------------------------------------------------

cat("\n8. VALIDAZIONE GLOBALE\n")
cat("------------------------------------------------------------\n")

checks <- data.frame(
  check = c(
    "No duplicate strict Ensembl keys",
    "Strict mapping symbols unique",
    "Strict mapping Ensembl unique",
    "All strict records have Ensembl ID",
    "All strict records have significant DE",
    "All strict records have consistent direction",
    "Unresolved records preserved",
    "Three previously missing genes recovered",
    "All records accounted for"
  ),
  passed = c(
    nrow(duplicate_rows) == 0,
    symbol_unique,
    ensembl_unique,
    all(
      !is.na(strict_rows$ensembl_id) &
        nzchar(strict_rows$ensembl_id)
    ),
    all(strict_rows$DE_significant_check),
    all(strict_rows$direction_consistency_check),
    nrow(unresolved_rows) ==
      sum(integrated$record_class == "UNRESOLVED"),
    nrow(recovery_audit) == 3 &&
      all(recovery_audit$recovery_validated),
    nrow(strict_rows) + nrow(unresolved_rows) ==
      nrow(integrated)
  ),
  stringsAsFactors = FALSE
)

print(checks)

write.csv(
  checks,
  file.path(output_dir, "global_validation_checks.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n9. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "all_core_gene_records_classified.csv",
  "strict_integrated_records.csv",
  "unresolved_records.csv",
  "duplicate_ensembl_core_gene_rows.csv",
  "strict_direction_consistency_audit.csv",
  "unresolved_mapping_audit.csv",
  "strict_symbol_uniqueness.csv",
  "strict_ensembl_uniqueness.csv",
  "three_gene_recovery_audit.csv",
  "core_gene_audit_summary.csv",
  "global_validation_checks.csv"
)

validation <- data.frame(
  file = output_files,
  exists = file.exists(
    file.path(output_dir, output_files)
  ),
  size_kb = round(
    file.info(
      file.path(output_dir, output_files)
    )$size / 1024,
    2
  )
)

print(validation)

if (!all(validation$exists)) {
  stop("Uno o più file di output non sono stati creati.")
}

# ------------------------------------------------------------
# 13. CONCLUSIONE
# ------------------------------------------------------------

cat("\n============================================================\n")

if (all(checks$passed)) {
  
  cat("CORE GENE AUDIT V3 COMPLETED\n")
  cat("GLOBAL VALIDATION PASSED\n")
  
} else {
  
  cat("CORE GENE AUDIT V3 COMPLETED WITH FAILED CHECKS\n")
  print(checks[!checks$passed, ])
  
}

cat("============================================================\n")

cat(
  "\nCartella output:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

if (!all(checks$passed)) {
  stop(
    "Validazione non superata. Consultare global_validation_checks.csv."
  )
}
