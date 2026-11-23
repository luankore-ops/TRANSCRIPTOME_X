
# ============================================================
# TRANSCRIPTOME-X - CORE GENE AUDIT
# Script: 19_core_gene_audit.R
#
# Obiettivi:
# 1. Verificare duplicazioni gene/comparison/direction
# 2. Verificare associazioni simbolo - Ensembl ID
# 3. Individuare geni senza corrispondenza DE significativa
# 4. Verificare coerenza della direzione
# 5. Rivalutare le classifiche dei core genes
#
# Non modifica gli output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE AUDIT\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

input_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_analysis"
)

output_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  "top_20_core_genes.csv",
  "core_gene_frequency_overall.csv",
  "core_genes_with_differential_expression.csv",
  "ambiguous_gene_symbols.csv",
  "core_gene_summary.csv"
)

missing_files <- required_files[
  !file.exists(file.path(input_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "File mancanti: ",
    paste(missing_files, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE DEI FILE\n")
cat("------------------------------------------------------------\n")

top20 <- read.csv(
  file.path(input_dir, "top_20_core_genes.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

frequency <- read.csv(
  file.path(input_dir, "core_gene_frequency_overall.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

integrated <- read.csv(
  file.path(input_dir, "core_genes_with_differential_expression.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ambiguous <- read.csv(
  file.path(input_dir, "ambiguous_gene_symbols.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_original <- read.csv(
  file.path(input_dir, "core_gene_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Top 20:", nrow(top20), "righe\n")
cat("Frequency:", nrow(frequency), "righe\n")
cat("Integrated:", nrow(integrated), "righe\n")
cat("Ambiguous symbols:", nrow(ambiguous), "righe\n")

# Verifica colonne necessarie

required_columns <- c(
  "comparison", "direction", "gene_name",
  "gene_id", "logFC", "FDR",
  "de_significant", "de_direction",
  "direction_consistent",
  "total_term_associations",
  "unique_terms",
  "best_adjusted_p_value"
)

missing_columns <- setdiff(
  required_columns,
  names(integrated)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti nel file integrato: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. CONTROLLO DUPLICATI
# ------------------------------------------------------------

cat("\n2. CONTROLLO DUPLICATI\n")
cat("------------------------------------------------------------\n")

integrated$key <- paste(
  integrated$comparison,
  integrated$direction,
  integrated$gene_name,
  sep = "|||"
)

duplicate_rows <- integrated[
  duplicated(integrated$key) |
    duplicated(integrated$key, fromLast = TRUE),
]

cat(
  "Righe coinvolte in duplicati comparison/direction/gene:",
  nrow(duplicate_rows),
  "\n"
)

write.csv(
  duplicate_rows,
  file.path(output_dir, "duplicate_core_gene_rows.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. AUDIT SIMBOLO - GENE ID
# ------------------------------------------------------------

cat("\n3. AUDIT SIMBOLO - GENE ID\n")
cat("------------------------------------------------------------\n")

group_key <- paste(
  integrated$comparison,
  integrated$direction,
  integrated$gene_name,
  sep = "|||"
)

groups <- split(
  seq_len(nrow(integrated)),
  group_key
)

symbol_audit <- do.call(
  rbind,
  lapply(groups, function(idx) {
    
    x <- integrated[idx, , drop = FALSE]
    
    ids <- unique(
      x$gene_id[
        !is.na(x$gene_id) &
          nzchar(x$gene_id)
      ]
    )
    
    data.frame(
      comparison = x$comparison[1],
      direction = x$direction[1],
      gene_name = x$gene_name[1],
      gene_id_count = length(ids),
      gene_ids = paste(ids, collapse = ";"),
      row_count = nrow(x),
      stringsAsFactors = FALSE
    )
  })
)

rownames(symbol_audit) <- NULL

symbol_audit$ambiguous <- (
  symbol_audit$gene_id_count > 1
)

cat(
  "Gruppi gene/comparison/direction:",
  nrow(symbol_audit),
  "\n"
)

cat(
  "Gruppi con più Gene ID:",
  sum(symbol_audit$ambiguous),
  "\n"
)

write.csv(
  symbol_audit,
  file.path(output_dir, "symbol_gene_id_audit.csv"),
  row.names = FALSE
)

# Confronto con il precedente report delle ambiguità

ambiguous$key <- paste(
  ambiguous$comparison,
  ambiguous$gene_name,
  sep = "|||"
)

symbol_audit$ambiguity_key <- paste(
  symbol_audit$comparison,
  symbol_audit$gene_name,
  sep = "|||"
)

symbol_audit$reported_as_ambiguous <- (
  symbol_audit$ambiguity_key %in% ambiguous$key
)

ambiguity_comparison <- aggregate(
  cbind(
    ambiguous = as.integer(symbol_audit$ambiguous),
    reported = as.integer(symbol_audit$reported_as_ambiguous)
  ) ~ comparison,
  data = symbol_audit,
  FUN = sum
)

write.csv(
  ambiguity_comparison,
  file.path(output_dir, "ambiguity_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. GENI SENZA DE SIGNIFICATIVA
# ------------------------------------------------------------

cat("\n4. GENI SENZA DE SIGNIFICATIVA\n")
cat("------------------------------------------------------------\n")

not_significant <- integrated[
  is.na(integrated$de_significant) |
    integrated$de_significant == FALSE,
]

cat(
  "Geni non classificati come DE significativi:",
  nrow(not_significant),
  "\n"
)

write.csv(
  not_significant,
  file.path(output_dir, "core_genes_without_significant_DE.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. CONTROLLO DIREZIONE
# ------------------------------------------------------------

cat("\n5. CONTROLLO COERENZA DELLA DIREZIONE\n")
cat("------------------------------------------------------------\n")

direction_audit <- integrated

direction_audit$expected_logFC_direction <- ifelse(
  direction_audit$logFC > 0,
  "UP",
  ifelse(direction_audit$logFC < 0, "DOWN", "ZERO")
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
  direction_audit$de_significant == TRUE &
    direction_audit$DE_direction_matches
)

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
  "Direzione coerente secondo controllo indipendente:",
  sum(
    direction_audit$direction_consistent_check,
    na.rm = TRUE
  ),
  "\n"
)

write.csv(
  direction_audit,
  file.path(output_dir, "direction_consistency_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. RIEPILOGO PER CONFRONTO E DIREZIONE
# ------------------------------------------------------------

cat("\n6. RIEPILOGO AUDIT\n")
cat("------------------------------------------------------------\n")

audit_summary <- do.call(
  rbind,
  lapply(
    split(
      seq_len(nrow(integrated)),
      paste(
        integrated$comparison,
        integrated$direction,
        sep = "|||"
      )
    ),
    function(idx) {
      
      x <- integrated[idx, , drop = FALSE]
      
      data.frame(
        comparison = x$comparison[1],
        direction = x$direction[1],
        integrated_rows = nrow(x),
        distinct_symbols = length(unique(x$gene_name)),
        distinct_gene_ids = length(unique(
          x$gene_id[
            !is.na(x$gene_id) &
              nzchar(x$gene_id)
          ]
        )),
        significant_DE = sum(
          x$de_significant == TRUE,
          na.rm = TRUE
        ),
        direction_consistent = sum(
          x$direction_consistent == TRUE,
          na.rm = TRUE
        ),
        ambiguous_symbols = sum(
          symbol_audit$ambiguous[
            symbol_audit$comparison == x$comparison[1] &
              symbol_audit$direction == x$direction[1]
          ],
          na.rm = TRUE
        )
      )
    }
  )
)

rownames(audit_summary) <- NULL

print(audit_summary)

write.csv(
  audit_summary,
  file.path(output_dir, "core_gene_audit_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. CLASSIFICA DEI CORE GENES NON AMBIGUI
# ------------------------------------------------------------

cat("\n7. CLASSIFICA CORE GENES NON AMBIGUI\n")
cat("------------------------------------------------------------\n")

unambiguous_keys <- symbol_audit[
  !symbol_audit$ambiguous,
  c("comparison", "direction", "gene_name")
]

unambiguous <- merge(
  integrated,
  unambiguous_keys,
  by = c("comparison", "direction", "gene_name"),
  all = FALSE,
  sort = FALSE
)

unambiguous <- unambiguous[
  order(
    unambiguous$comparison,
    unambiguous$direction,
    -unambiguous$unique_terms,
    unambiguous$best_adjusted_p_value
  ),
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
]

write.csv(
  top_unambiguous,
  file.path(output_dir, "top_20_unambiguous_core_genes.csv"),
  row.names = FALSE
)

cat(
  "Geni non ambigui:",
  nrow(unambiguous),
  "\n"
)

cat(
  "Top genes non ambigui:",
  nrow(top_unambiguous),
  "\n"
)

# ------------------------------------------------------------
# 9. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n8. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "duplicate_core_gene_rows.csv",
  "symbol_gene_id_audit.csv",
  "ambiguity_summary.csv",
  "core_genes_without_significant_DE.csv",
  "direction_consistency_audit.csv",
  "core_gene_audit_summary.csv",
  "unambiguous_core_genes_ranked.csv",
  "top_20_unambiguous_core_genes.csv"
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
  stop("Validazione fallita: uno o più file non sono stati creati.")
}

# ------------------------------------------------------------
# 10. CONCLUSIONE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE AUDIT COMPLETED\n")
cat("============================================================\n")

cat(
  "\nCartella risultati:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nCORE GENE AUDIT VALIDATED.\n")
