
# ============================================================
# TRANSCRIPTOME-X - CORE GENE AUDIT V2
# Script: 21_core_gene_audit_ensembl.R
#
# Audit basato su identificativi Ensembl stabili.
#
# Obiettivi:
# 1. Verificare duplicazioni gene/comparison/direction
# 2. Verificare unicita e coerenza degli identificativi Ensembl
# 3. Verificare mapping strict e simboli ambigui
# 4. Verificare integrazione con differential expression
# 5. Verificare significativita DE e coerenza della direzione
# 6. Verificare il recupero dei tre geni precedentemente persi
# 7. Rivalutare le classifiche dei core genes
#
# Non modifica gli script e gli output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE AUDIT V2\n")
cat("ENSEMBL-BASED VALIDATION\n")
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
  "reports", "core_gene_audit_v2"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  "gene_pathway_membership_ensembl.csv",
  "gene_symbol_mapping_audit.csv",
  "gene_symbol_mapping_strict.csv",
  "core_gene_frequency_overall.csv",
  "core_genes_with_differential_expression.csv",
  "top_20_core_genes.csv",
  "core_gene_summary.csv"
)

missing_files <- required_files[
  !file.exists(file.path(input_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "File mancanti nella cartella V2: ",
    paste(missing_files, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. FUNZIONI DI SUPPORTO
# ------------------------------------------------------------

read_project_csv <- function(filename) {
  read.csv(
    file.path(input_dir, filename),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

check_columns <- function(data, required, object_name) {
  
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      "Colonne mancanti in ", object_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
}

is_true <- function(x) {
  !is.na(x) & as.logical(x)
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\.[0-9]+$", "", x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

safe_unique_count <- function(x) {
  length(unique(x[!is.na(x) & nzchar(as.character(x))]))
}

# ------------------------------------------------------------
# 3. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE DEI FILE\n")
cat("------------------------------------------------------------\n")

membership <- read_project_csv(
  "gene_pathway_membership_ensembl.csv"
)

mapping_audit <- read_project_csv(
  "gene_symbol_mapping_audit.csv"
)

mapping_strict <- read_project_csv(
  "gene_symbol_mapping_strict.csv"
)

frequency <- read_project_csv(
  "core_gene_frequency_overall.csv"
)

integrated <- read_project_csv(
  "core_genes_with_differential_expression.csv"
)

top20_original <- read_project_csv(
  "top_20_core_genes.csv"
)

summary_original <- read_project_csv(
  "core_gene_summary.csv"
)

cat("Gene membership:", nrow(membership), "righe\n")
cat("Mapping audit:", nrow(mapping_audit), "righe\n")
cat("Mapping strict:", nrow(mapping_strict), "righe\n")
cat("Frequency:", nrow(frequency), "righe\n")
cat("Integrated:", nrow(integrated), "righe\n")
cat("Top 20 originale:", nrow(top20_original), "righe\n")

# ------------------------------------------------------------
# 4. VALIDAZIONE DELLE COLONNE
# ------------------------------------------------------------

cat("\n2. VALIDAZIONE DELLE COLONNE\n")
cat("------------------------------------------------------------\n")

required_integrated <- c(
  "comparison",
  "direction",
  "enrichment_symbol",
  "ensembl_id",
  "mapping_status",
  "de_gene_name",
  "logFC",
  "FDR",
  "de_significant",
  "de_direction",
  "direction_consistent",
  "total_term_associations",
  "unique_terms",
  "best_adjusted_p_value"
)

check_columns(
  integrated,
  required_integrated,
  "core_genes_with_differential_expression.csv"
)

check_columns(
  membership,
  c(
    "comparison",
    "direction",
    "enrichment_symbol",
    "ensembl_id",
    "mapping_status",
    "gene_key"
  ),
  "gene_pathway_membership_ensembl.csv"
)

check_columns(
  mapping_audit,
  c("SYMBOL", "mapping_status"),
  "gene_symbol_mapping_audit.csv"
)

check_columns(
  mapping_strict,
  c("SYMBOL", "ENSEMBL"),
  "gene_symbol_mapping_strict.csv"
)

# Normalizzazione degli identificativi Ensembl.
# Elimina eventuali suffissi di versione.

integrated$ensembl_id <- clean_id(
  integrated$ensembl_id
)

membership$ensembl_id <- clean_id(
  membership$ensembl_id
)

mapping_strict$ENSEMBL <- clean_id(
  mapping_strict$ENSEMBL
)

cat("Colonne obbligatorie: OK\n")
cat("Identificativi Ensembl normalizzati.\n")

# ------------------------------------------------------------
# 5. CONTROLLO DEI DUPLICATI
# ------------------------------------------------------------

cat("\n3. CONTROLLO DUPLICATI\n")
cat("------------------------------------------------------------\n")

integrated$key <- paste(
  integrated$comparison,
  integrated$direction,
  integrated$ensembl_id,
  sep = "|||"
)

integrated$key[
  is.na(integrated$ensembl_id)
] <- paste0(
  "UNMAPPED|||",
  which(is.na(integrated$ensembl_id))
)

duplicate_rows <- integrated[
  duplicated(integrated$key) |
    duplicated(integrated$key, fromLast = TRUE),
  ,
  drop = FALSE
]

cat(
  "Righe coinvolte in duplicati comparison/direction/Ensembl:",
  nrow(duplicate_rows),
  "\n"
)

write.csv(
  duplicate_rows,
  file.path(output_dir, "duplicate_ensembl_core_gene_rows.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. AUDIT ENSEMBL E SIMBOLI
# ------------------------------------------------------------

cat("\n4. AUDIT ENSEMBL E SIMBOLI\n")
cat("------------------------------------------------------------\n")

valid_ensembl <- integrated[
  !is.na(integrated$ensembl_id),
  ,
  drop = FALSE
]

ensembl_groups <- split(
  seq_len(nrow(valid_ensembl)),
  paste(
    valid_ensembl$comparison,
    valid_ensembl$direction,
    valid_ensembl$ensembl_id,
    sep = "|||"
  )
)

ensembl_audit <- do.call(
  rbind,
  lapply(ensembl_groups, function(idx) {
    
    x <- valid_ensembl[idx, , drop = FALSE]
    
    data.frame(
      comparison = x$comparison[1],
      direction = x$direction[1],
      ensembl_id = x$ensembl_id[1],
      enrichment_symbols = paste(
        sort(unique(na.omit(x$enrichment_symbol))),
        collapse = ";"
      ),
      de_gene_names = paste(
        sort(unique(na.omit(x$de_gene_name))),
        collapse = ";"
      ),
      enrichment_symbol_count = safe_unique_count(
        x$enrichment_symbol
      ),
      de_gene_name_count = safe_unique_count(
        x$de_gene_name
      ),
      row_count = nrow(x),
      stringsAsFactors = FALSE
    )
  })
)

rownames(ensembl_audit) <- NULL

ensembl_audit$multiple_enrichment_symbols <- (
  ensembl_audit$enrichment_symbol_count > 1
)

ensembl_audit$multiple_de_symbols <- (
  ensembl_audit$de_gene_name_count > 1
)

cat(
  "Gruppi Ensembl/comparison/direction:",
  nrow(ensembl_audit),
  "\n"
)

cat(
  "Gruppi con più simboli enrichment:",
  sum(ensembl_audit$multiple_enrichment_symbols),
  "\n"
)

cat(
  "Gruppi con più simboli DE:",
  sum(ensembl_audit$multiple_de_symbols),
  "\n"
)

write.csv(
  ensembl_audit,
  file.path(output_dir, "ensembl_symbol_consistency_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. AUDIT MAPPING STRICT E AMBIGUITA
# ------------------------------------------------------------

cat("\n5. AUDIT MAPPING STRICT E AMBIGUITA\n")
cat("------------------------------------------------------------\n")

mapping_status_summary <- as.data.frame(
  table(mapping_audit$mapping_status),
  stringsAsFactors = FALSE
)

names(mapping_status_summary) <- c(
  "mapping_status",
  "symbol_count"
)

print(mapping_status_summary)

write.csv(
  mapping_status_summary,
  file.path(output_dir, "mapping_status_audit_summary.csv"),
  row.names = FALSE
)

# Simboli ambigui: non devono essere integrati
# attraverso un mapping Ensembl arbitrario.

ambiguous_symbols <- mapping_audit$SYMBOL[
  mapping_audit$mapping_status == "AMBIGUOUS_SYMBOL"
]

ambiguous_symbols <- unique(
  ambiguous_symbols[!is.na(ambiguous_symbols)]
)

ambiguous_in_integrated <- integrated[
  integrated$enrichment_symbol %in% ambiguous_symbols,
  ,
  drop = FALSE
]

cat(
  "Simboli ambigui nel mapping:",
  length(ambiguous_symbols),
  "\n"
)

cat(
  "Righe integrate con simbolo ambiguo:",
  nrow(ambiguous_in_integrated),
  "\n"
)

write.csv(
  ambiguous_in_integrated,
  file.path(output_dir, "ambiguous_symbols_integration_audit.csv"),
  row.names = FALSE
)

# Verifica che i mapping strict siano realmente univoci.

strict_mapping <- mapping_strict[
  !is.na(mapping_strict$SYMBOL) &
    !is.na(mapping_strict$ENSEMBL),
  ,
  drop = FALSE
]

strict_symbol_counts <- aggregate(
  ENSEMBL ~ SYMBOL,
  data = strict_mapping,
  FUN = function(x) length(unique(x))
)

strict_ensembl_counts <- aggregate(
  SYMBOL ~ ENSEMBL,
  data = strict_mapping,
  FUN = function(x) length(unique(x))
)

strict_symbol_counts$unique_ensembl <- (
  strict_symbol_counts$ENSEMBL == 1
)

strict_ensembl_counts$unique_symbol <- (
  strict_ensembl_counts$SYMBOL == 1
)

cat(
  "Simboli strict con un solo Ensembl:",
  sum(strict_symbol_counts$unique_ensembl),
  "\n"
)

cat(
  "Ensembl strict associati a un solo simbolo:",
  sum(strict_ensembl_counts$unique_symbol),
  "\n"
)

write.csv(
  strict_symbol_counts,
  file.path(output_dir, "strict_symbol_uniqueness_audit.csv"),
  row.names = FALSE
)

write.csv(
  strict_ensembl_counts,
  file.path(output_dir, "strict_ensembl_uniqueness_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. GENI SENZA DE SIGNIFICATIVA
# ------------------------------------------------------------

cat("\n6. AUDIT DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")

integrated$de_significant_check <- is_true(
  integrated$de_significant
)

not_significant <- integrated[
  !integrated$de_significant_check,
  ,
  drop = FALSE
]

cat(
  "Righe senza DE significativa:",
  nrow(not_significant),
  "\n"
)

write.csv(
  not_significant,
  file.path(output_dir, "core_genes_without_significant_DE.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. CONTROLLO COERENZA DELLA DIREZIONE
# ------------------------------------------------------------

cat("\n7. CONTROLLO COERENZA DELLA DIREZIONE\n")
cat("------------------------------------------------------------\n")

direction_audit <- integrated

direction_audit$expected_logFC_direction <- ifelse(
  is.na(direction_audit$logFC),
  NA_character_,
  ifelse(
    direction_audit$logFC > 0,
    "UP",
    ifelse(direction_audit$logFC < 0, "DOWN", "ZERO")
  )
)

direction_audit$logFC_matches_direction <- (
  direction_audit$expected_logFC_direction ==
    direction_audit$direction
)

direction_audit$DE_direction_matches <- (
  direction_audit$de_direction ==
    direction_audit$direction
)

direction_audit$direction_consistent_check <- (
  direction_audit$de_significant_check &
    direction_audit$logFC_matches_direction &
    direction_audit$DE_direction_matches
)

direction_audit$direction_consistent_check[
  is.na(direction_audit$direction_consistent_check)
] <- FALSE

cat(
  "Direzione logFC non corrispondente:",
  sum(
    !direction_audit$logFC_matches_direction,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Direzione DE non corrispondente:",
  sum(
    !direction_audit$DE_direction_matches,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Righe con significativita e direzione coerenti:",
  sum(direction_audit$direction_consistent_check),
  "\n"
)

write.csv(
  direction_audit,
  file.path(output_dir, "direction_consistency_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. RIEPILOGO PER CONFRONTO E DIREZIONE
# ------------------------------------------------------------

cat("\n8. RIEPILOGO AUDIT\n")
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
    
    direction_x <- direction_audit[idx, , drop = FALSE]
    
    data.frame(
      comparison = x$comparison[1],
      direction = x$direction[1],
      integrated_rows = nrow(x),
      distinct_enrichment_symbols = safe_unique_count(
        x$enrichment_symbol
      ),
      distinct_ensembl_ids = safe_unique_count(
        x$ensembl_id
      ),
      matched_to_ensembl = sum(
        !is.na(x$ensembl_id)
      ),
      significant_DE = sum(
        x$de_significant_check,
        na.rm = TRUE
      ),
      direction_consistent = sum(
        direction_x$direction_consistent_check,
        na.rm = TRUE
      ),
      missing_ensembl = sum(
        is.na(x$ensembl_id)
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
# 11. CLASSIFICA DEI CORE GENES NON AMBIGUI
# ------------------------------------------------------------

cat("\n9. CLASSIFICA CORE GENES NON AMBIGUI\n")
cat("------------------------------------------------------------\n")

unambiguous <- integrated[
  !is.na(integrated$ensembl_id) &
    integrated$mapping_status == "STRICT_ONE_TO_ONE",
  ,
  drop = FALSE
]

# Esclude eventuali righe duplicate sulla stessa chiave.
unambiguous <- unambiguous[
  !duplicated(
    paste(
      unambiguous$comparison,
      unambiguous$direction,
      unambiguous$ensembl_id,
      sep = "|||"
    )
  ),
  ,
  drop = FALSE
]

unambiguous <- unambiguous[
  order(
    unambiguous$comparison,
    unambiguous$direction,
    -unambiguous$unique_terms,
    unambiguous$best_adjusted_p_value,
    -unambiguous$total_term_associations
  ),
  ,
  drop = FALSE
]

unambiguous$rank <- ave(
  seq_len(nrow(unambiguous)),
  interaction(
    unambiguous$comparison,
    unambiguous$direction,
    drop = TRUE
  ),
  FUN = seq_along
)

write.csv(
  unambiguous,
  file.path(output_dir, "unambiguous_core_genes_ranked.csv"),
  row.names = FALSE
)

top_unambiguous <- unambiguous[
  unambiguous$rank <= 20,
  ,
  drop = FALSE
]

write.csv(
  top_unambiguous,
  file.path(output_dir, "top_20_unambiguous_core_genes.csv"),
  row.names = FALSE
)

cat(
  "Core genes strict integrati:",
  nrow(unambiguous),
  "\n"
)

cat(
  "Righe nelle classifiche Top 20:",
  nrow(top_unambiguous),
  "\n"
)

# ------------------------------------------------------------
# 12. AUDIT RECUPERO DEI TRE GENI
# ------------------------------------------------------------

cat("\n10. AUDIT RECUPERO DEI TRE GENI\n")
cat("------------------------------------------------------------\n")

target_genes <- data.frame(
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
  expected_ensembl_id = c(
    "ENSG00000188522",
    "ENSG00000112697",
    "ENSG00000146540"
  ),
  stringsAsFactors = FALSE
)

target_audit <- merge(
  target_genes,
  integrated,
  by = c("comparison", "enrichment_symbol"),
  all.x = TRUE,
  sort = FALSE
)

target_audit$ensembl_match <- (
  target_audit$ensembl_id ==
    target_audit$expected_ensembl_id
)

target_audit$de_symbol_match <- (
  target_audit$de_gene_name ==
    target_audit$expected_de_gene_name
)

target_audit$significant_DE_check <- is_true(
  target_audit$de_significant
)

target_audit$recovery_validated <- (
  target_audit$ensembl_match &
    target_audit$de_symbol_match &
    target_audit$significant_DE_check &
    target_audit$mapping_status == "STRICT_ONE_TO_ONE"
)

print(
  target_audit[
    ,
    intersect(
      c(
        "comparison",
        "enrichment_symbol",
        "de_gene_name",
        "expected_de_gene_name",
        "ensembl_id",
        "logFC",
        "FDR",
        "de_significant",
        "mapping_status",
        "recovery_validated"
      ),
      names(target_audit)
    ),
    drop = FALSE
  ]
)

write.csv(
  target_audit,
  file.path(output_dir, "three_gene_recovery_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. VALIDAZIONE GLOBALE
# ------------------------------------------------------------

cat("\n11. VALIDAZIONE GLOBALE\n")
cat("------------------------------------------------------------\n")

checks <- data.frame(
  check = c(
    "No duplicate Ensembl keys",
    "No ambiguous symbols integrated",
    "All integrated rows have Ensembl ID",
    "All integrated rows have significant DE",
    "All integrated rows have consistent direction",
    "Three previously missing genes recovered"
  ),
  passed = c(
    nrow(duplicate_rows) == 0,
    nrow(ambiguous_in_integrated) == 0,
    all(!is.na(integrated$ensembl_id)),
    all(integrated$de_significant_check),
    all(direction_audit$direction_consistent_check),
    nrow(target_audit) == 3 &&
      all(target_audit$recovery_validated)
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
# 14. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n12. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "duplicate_ensembl_core_gene_rows.csv",
  "ensembl_symbol_consistency_audit.csv",
  "mapping_status_audit_summary.csv",
  "ambiguous_symbols_integration_audit.csv",
  "strict_symbol_uniqueness_audit.csv",
  "strict_ensembl_uniqueness_audit.csv",
  "core_genes_without_significant_DE.csv",
  "direction_consistency_audit.csv",
  "core_gene_audit_summary.csv",
  "unambiguous_core_genes_ranked.csv",
  "top_20_unambiguous_core_genes.csv",
  "three_gene_recovery_audit.csv",
  "global_validation_checks.csv"
)

output_paths <- file.path(
  output_dir,
  output_files
)

validation <- data.frame(
  file = output_files,
  exists = file.exists(output_paths),
  size_kb = round(
    file.info(output_paths)$size / 1024,
    2
  )
)

print(validation)

if (!all(validation$exists)) {
  stop(
    "Validazione fallita: uno o più file non sono stati creati."
  )
}

if (!all(checks$passed)) {
  stop(
    "Audit completato, ma uno o più controlli globali sono falliti. ",
    "Consultare global_validation_checks.csv."
  )
}

# ------------------------------------------------------------
# 15. CONCLUSIONE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE AUDIT V2 COMPLETED\n")
cat("============================================================\n")

cat(
  "\nGeni integrati:",
  nrow(integrated),
  "\n"
)

cat(
  "Mapping Ensembl strict:",
  sum(!is.na(integrated$ensembl_id)),
  "\n"
)

cat(
  "Geni DE significativi:",
  sum(integrated$de_significant_check),
  "\n"
)

cat(
  "Geni con direzione coerente:",
  sum(direction_audit$direction_consistent_check),
  "\n"
)

cat(
  "Tre geni recuperati:",
  sum(target_audit$recovery_validated),
  "/ 3\n"
)

cat(
  "\nCartella risultati:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nCORE GENE AUDIT V2 VALIDATED.\n")
