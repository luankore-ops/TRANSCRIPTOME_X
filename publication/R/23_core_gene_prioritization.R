
# ============================================================
# TRANSCRIPTOME-X - CORE GENE PRIORITIZATION
# Script: 23_core_gene_prioritization.R
#
# Obiettivi:
# 1. Correggere il riepilogo della direzione di R22
# 2. Identificare i core genes ricorrenti nei due knockdown
# 3. Integrare pathway frequency e differential expression
# 4. Calcolare uno score composito di prioritizzazione
# 5. Generare classifiche e report
#
# Input: R22 - strict_integrated_records.csv
# Non modifica gli output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE PRIORITIZATION\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

input_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_audit_v3"
)

output_dir <- file.path(
  "data", "processed",
  "functional_enrichment_strict",
  "reports", "core_gene_prioritization"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  input_dir,
  "strict_integrated_records.csv"
)

if (!file.exists(input_file)) {
  stop("File di input non trovato: ", input_file)
}

# ------------------------------------------------------------
# 2. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE RECORD STRICT\n")
cat("------------------------------------------------------------\n")

genes <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Record importati:", nrow(genes), "\n")

required_columns <- c(
  "comparison",
  "direction",
  "enrichment_symbol",
  "ensembl_id",
  "de_gene_name",
  "logFC",
  "FDR",
  "de_significant",
  "de_direction",
  "total_term_associations",
  "unique_terms",
  "best_adjusted_p_value"
)

missing_columns <- setdiff(
  required_columns,
  names(genes)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. NORMALIZZAZIONE E CONTROLLO DATI
# ------------------------------------------------------------

cat("\n2. NORMALIZZAZIONE E CONTROLLO\n")
cat("------------------------------------------------------------\n")

genes$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  trimws(as.character(genes$ensembl_id))
)

numeric_columns <- c(
  "logFC",
  "FDR",
  "total_term_associations",
  "unique_terms",
  "best_adjusted_p_value"
)

for (column in numeric_columns) {
  genes[[column]] <- suppressWarnings(
    as.numeric(genes[[column]])
  )
}

genes$de_significant <- as.logical(
  genes$de_significant
)

if (anyNA(genes$ensembl_id) ||
    any(!nzchar(genes$ensembl_id))) {
  stop("Sono presenti identificativi Ensembl mancanti.")
}

if (anyDuplicated(
  genes[c("comparison", "direction", "ensembl_id")]
) > 0) {
  stop("Sono presenti duplicati nella chiave Ensembl.")
}

if (!all(genes$de_significant %in% TRUE)) {
  stop("Sono presenti record senza DE significativa.")
}

# Ricostruzione indipendente della coerenza direzionale

genes$expected_direction <- ifelse(
  genes$logFC > 0,
  "UP",
  ifelse(genes$logFC < 0, "DOWN", "ZERO")
)

genes$logFC_direction_match <- (
  genes$expected_direction == genes$direction
)

genes$DE_direction_match <- (
  genes$de_direction == genes$direction
)

genes$direction_consistent_check <- (
  genes$de_significant &
    genes$logFC_direction_match &
    genes$DE_direction_match
)

cat(
  "Record con direzione coerente:",
  sum(genes$direction_consistent_check, na.rm = TRUE),
  "/",
  nrow(genes),
  "\n"
)

if (!all(genes$direction_consistent_check)) {
  stop("La verifica indipendente della direzione non è superata.")
}

# ------------------------------------------------------------
# 4. RIEPILOGO R22 CORRETTO
# ------------------------------------------------------------

cat("\n3. RIEPILOGO AUDIT CORRETTO\n")
cat("------------------------------------------------------------\n")

summary_groups <- split(
  seq_len(nrow(genes)),
  paste(
    genes$comparison,
    genes$direction,
    sep = "|||"
  )
)

corrected_audit_summary <- do.call(
  rbind,
  lapply(summary_groups, function(idx) {
    
    x <- genes[idx, , drop = FALSE]
    
    data.frame(
      comparison = x$comparison[1],
      direction = x$direction[1],
      strict_records = nrow(x),
      distinct_ensembl = length(unique(x$ensembl_id)),
      significant_DE = sum(x$de_significant, na.rm = TRUE),
      direction_consistent = sum(
        x$direction_consistent_check,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
)

rownames(corrected_audit_summary) <- NULL

print(corrected_audit_summary)

write.csv(
  corrected_audit_summary,
  file.path(output_dir, "corrected_core_gene_audit_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. FUNZIONE PER SCORE PERCENTILE
# ------------------------------------------------------------

percentile_score <- function(x) {
  
  x <- as.numeric(x)
  
  result <- rep(NA_real_, length(x))
  
  valid <- which(is.finite(x))
  
  if (length(valid) == 0) {
    return(result)
  }
  
  if (length(valid) == 1) {
    result[valid] <- 0.5
    return(result)
  }
  
  result[valid] <- (
    rank(x[valid], ties.method = "average") - 1
  ) / (length(valid) - 1)
  
  result
}

# ------------------------------------------------------------
# 6. AGGREGAZIONE PER GENE ENSEMBL
# ------------------------------------------------------------

cat("\n4. AGGREGAZIONE PER GENE\n")
cat("------------------------------------------------------------\n")

gene_groups <- split(
  seq_len(nrow(genes)),
  genes$ensembl_id
)

gene_summary <- do.call(
  rbind,
  lapply(gene_groups, function(idx) {
    
    x <- genes[idx, , drop = FALSE]
    
    comparisons <- unique(x$comparison)
    directions <- unique(x$direction)
    
    comparison_count <- length(comparisons)
    
    cross_status <- if (comparison_count == 2 &&
                        length(directions) == 1) {
      
      "BOTH_KNOCKDOWNS_SAME_DIRECTION"
      
    } else if (comparison_count == 2 &&
               length(directions) > 1) {
      
      "OPPOSITE_DIRECTIONS"
      
    } else if (comparison_count == 1) {
      
      "SINGLE_KNOCKDOWN"
      
    } else {
      
      "OTHER"
    }
    
    data.frame(
      ensembl_id = x$ensembl_id[1],
      
      enrichment_symbols = paste(
        sort(unique(x$enrichment_symbol)),
        collapse = ";"
      ),
      
      de_gene_names = paste(
        sort(unique(x$de_gene_name)),
        collapse = ";"
      ),
      
      comparisons = paste(
        sort(comparisons),
        collapse = ";"
      ),
      
      directions = paste(
        sort(directions),
        collapse = ";"
      ),
      
      n_comparisons = comparison_count,
      
      cross_knockdown_status = cross_status,
      
      total_unique_terms_across_comparisons = sum(
        x$unique_terms,
        na.rm = TRUE
      ),
      
      total_term_associations = sum(
        x$total_term_associations,
        na.rm = TRUE
      ),
      
      best_DE_FDR = min(
        x$FDR,
        na.rm = TRUE
      ),
      
      best_abs_logFC = max(
        abs(x$logFC),
        na.rm = TRUE
      ),
      
      max_abs_logFC = max(
        abs(x$logFC),
        na.rm = TRUE
      ),
      
      all_DE_significant = all(
        x$de_significant
      ),
      
      all_directions_consistent = all(
        x$direction_consistent_check
      ),
      
      stringsAsFactors = FALSE
    )
  })
)

rownames(gene_summary) <- NULL

cat(
  "Geni Ensembl distinti:",
  nrow(gene_summary),
  "\n"
)

# ------------------------------------------------------------
# 7. CALCOLO SCORE DI PRIORITIZZAZIONE
# ------------------------------------------------------------

cat("\n5. CALCOLO PRIORITY SCORE\n")
cat("------------------------------------------------------------\n")

# 1. Ricorrenza nei due knockdown
gene_summary$recurrence_score <- (
  gene_summary$n_comparisons / 2
)

# 2. Ampiezza della presenza nei pathway
gene_summary$pathway_breadth_score <- percentile_score(
  gene_summary$total_unique_terms_across_comparisons
)

# 3. Significatività DE
gene_summary$DE_significance_score <- percentile_score(
  -log10(
    pmax(gene_summary$best_DE_FDR, 1e-300)
  )
)

# 4. Ampiezza della variazione di espressione
gene_summary$effect_size_score <- percentile_score(
  gene_summary$best_abs_logFC
)

# Score composito:
# 35% ricorrenza
# 25% pathway breadth
# 20% significatività DE
# 20% ampiezza logFC

gene_summary$priority_score <- (
  0.35 * gene_summary$recurrence_score +
    0.25 * gene_summary$pathway_breadth_score +
    0.20 * gene_summary$DE_significance_score +
    0.20 * gene_summary$effect_size_score
)

gene_summary <- gene_summary[
  order(
    -gene_summary$priority_score,
    -gene_summary$n_comparisons,
    gene_summary$best_DE_FDR
  ),
]

gene_summary$priority_rank <- seq_len(
  nrow(gene_summary)
)

gene_summary$priority_tier <- ifelse(
  gene_summary$cross_knockdown_status ==
    "BOTH_KNOCKDOWNS_SAME_DIRECTION",
  "HIGH_REPLICATION",
  ifelse(
    gene_summary$cross_knockdown_status ==
      "OPPOSITE_DIRECTIONS",
    "DIRECTION_DISCORDANT",
    "SINGLE_KNOCKDOWN"
  )
)

write.csv(
  gene_summary,
  file.path(output_dir, "core_gene_priority_ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. GENI CONCORDANTI TRA I DUE KNOCKDOWN
# ------------------------------------------------------------

cat("\n6. CORE GENES CONCORDANTI\n")
cat("------------------------------------------------------------\n")

concordant_genes <- gene_summary[
  gene_summary$cross_knockdown_status ==
    "BOTH_KNOCKDOWNS_SAME_DIRECTION",
  ,
  drop = FALSE
]

concordant_genes <- concordant_genes[
  order(
    -concordant_genes$priority_score
  ),
]

write.csv(
  concordant_genes,
  file.path(output_dir, "concordant_core_genes.csv"),
  row.names = FALSE
)

cat(
  "Geni presenti in entrambi i knockdown con direzione concordante:",
  nrow(concordant_genes),
  "\n"
)

# ------------------------------------------------------------
# 9. CLASSIFICHE PER CONFRONTO E DIREZIONE
# ------------------------------------------------------------

cat("\n7. CLASSIFICHE PER CONFRONTO E DIREZIONE\n")
cat("------------------------------------------------------------\n")

genes$DE_significance_score <- percentile_score(
  -log10(pmax(genes$FDR, 1e-300))
)

genes$pathway_breadth_score <- percentile_score(
  genes$unique_terms
)

genes$effect_size_score <- percentile_score(
  abs(genes$logFC)
)

genes$comparison_priority_score <- (
  0.40 * genes$pathway_breadth_score +
    0.30 * genes$DE_significance_score +
    0.30 * genes$effect_size_score
)

genes <- genes[
  order(
    genes$comparison,
    genes$direction,
    -genes$comparison_priority_score
  ),
]

ranking_groups <- split(
  seq_len(nrow(genes)),
  paste(
    genes$comparison,
    genes$direction,
    sep = "|||"
  )
)

ranked_records <- do.call(
  rbind,
  lapply(ranking_groups, function(idx) {
    
    x <- genes[idx, , drop = FALSE]
    
    x$comparison_rank <- rank(
      -x$comparison_priority_score,
      ties.method = "min"
    )
    
    x
  })
)

rownames(ranked_records) <- NULL

write.csv(
  ranked_records,
  file.path(output_dir, "comparison_direction_rankings.csv"),
  row.names = FALSE
)

top20 <- ranked_records[
  ranked_records$comparison_rank <= 20,
]

write.csv(
  top20,
  file.path(output_dir, "top_20_per_comparison_direction.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. REPORT DEI RISULTATI
# ------------------------------------------------------------

cat("\n8. REPORT\n")
cat("------------------------------------------------------------\n")

report_file <- file.path(
  output_dir,
  "TRANSCRIPTOME_X_CORE_GENE_PRIORITIZATION_REPORT.md"
)

report_lines <- c(
  "# TRANSCRIPTOME-X - Core Gene Prioritization",
  "",
  "## Metodo",
  "",
  "La prioritizzazione integra quattro componenti:",
  "",
  "- Ricorrenza nei due knockdown: 35%",
  "- Ampiezza dei pathway: 25%",
  "- Significatività DE: 20%",
  "- Ampiezza del logFC: 20%",
  "",
  "Gli score sono percentili relativi alla distribuzione dei geni analizzati.",
  "",
  "## Risultati",
  "",
  paste("- Record strict analizzati:", nrow(genes)),
  paste("- Geni Ensembl distinti:", nrow(gene_summary)),
  paste("- Geni concordanti tra i due knockdown:", nrow(concordant_genes)),
  paste("- Geni nelle classifiche Top 20:", nrow(top20)),
  "",
  "## Interpretazione",
  "",
  "Il priority score è uno strumento esplorativo di ordinamento.",
  "Non rappresenta una misura diretta di causalità o importanza biologica.",
  "La ricorrenza nei pathway non dimostra attivazione o inibizione funzionale.",
  "La concordanza tra knockdown indica consistenza osservata, non causalità.",
  "",
  "## Primi 20 geni prioritari",
  ""
)

top_display <- head(gene_summary, 20)

for (i in seq_len(nrow(top_display))) {
  
  report_lines <- c(
    report_lines,
    paste0(
      i, ". ",
      top_display$de_gene_names[i],
      " (", top_display$ensembl_id[i], ")",
      " - score: ",
      round(top_display$priority_score[i], 4),
      " - ",
      top_display$cross_knockdown_status[i]
    )
  )
}

writeLines(
  report_lines,
  report_file
)

# ------------------------------------------------------------
# 11. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n9. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "corrected_core_gene_audit_summary.csv",
  "core_gene_priority_ranking.csv",
  "concordant_core_genes.csv",
  "comparison_direction_rankings.csv",
  "top_20_per_comparison_direction.csv",
  "TRANSCRIPTOME_X_CORE_GENE_PRIORITIZATION_REPORT.md"
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
  stop("Validazione output fallita.")
}

# ------------------------------------------------------------
# 12. RIEPILOGO FINALE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORE GENE PRIORITIZATION COMPLETED\n")
cat("============================================================\n")

cat(
  "\nGeni prioritizzati:",
  nrow(gene_summary),
  "\n"
)

cat(
  "Geni concordanti tra i due knockdown:",
  nrow(concordant_genes),
  "\n"
)

cat("\nTOP 10 CORE GENES\n")
cat("------------------------------------------------------------\n")

print(
  head(
    gene_summary[
      ,
      c(
        "priority_rank",
        "de_gene_names",
        "n_comparisons",
        "cross_knockdown_status",
        "priority_score"
      )
    ],
    10
  ),
  row.names = FALSE
)

cat(
  "\nCartella output:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nCORE GENE PRIORITIZATION VALIDATED.\n")
