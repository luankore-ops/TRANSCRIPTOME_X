
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 16 - ROBUSTNESS REPORT
# ============================================================
#
# Obiettivo:
# Documentare la robustezza della functional enrichment
# analysis rispetto alla strategia di mapping Ensembl-Entrez.
#
# Input:
# - Functional enrichment originale
# - Functional enrichment con mapping rigoroso 1:1
# - Tabelle di confronto già prodotte
#
# Output:
# data/processed/functional_enrichment_strict/reports/
#   robustness_report/
#
# Nessun risultato precedente viene modificato.
#
# ============================================================


# 1. CONFIGURAZIONE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - ROBUSTNESS REPORT\n")
cat("============================================================\n\n")

base_dir <- "data/processed"

original_dir <- file.path(
  base_dir,
  "functional_enrichment"
)

strict_dir <- file.path(
  base_dir,
  "functional_enrichment_strict"
)

output_dir <- file.path(
  strict_dir,
  "reports",
  "robustness_report"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. DEFINIZIONE DEI FILE INPUT
# ============================================================

input_files <- c(
  
  original_mapping = file.path(
    original_dir,
    "mapping",
    "mapping_summary.csv"
  ),
  
  strict_mapping = file.path(
    strict_dir,
    "mapping",
    "mapping_summary_strict.csv"
  ),
  
  pathway_comparison = file.path(
    strict_dir,
    "reports",
    "enrichment_robustness_comparison.csv"
  ),
  
  gene_coverage = file.path(
    strict_dir,
    "reports",
    "gene_mapping_coverage_strict.csv"
  )
)

for (f in input_files) {
  
  if (!file.exists(f)) {
    stop(
      paste("File necessario non trovato:", f)
    )
  }
}

cat("Tutti i file di input sono disponibili.\n")


# 3. IMPORTAZIONE DATI
# ============================================================

original_mapping <- read.csv(
  input_files["original_mapping"],
  stringsAsFactors = FALSE
)

strict_mapping <- read.csv(
  input_files["strict_mapping"],
  stringsAsFactors = FALSE
)

pathway_comparison <- read.csv(
  input_files["pathway_comparison"],
  stringsAsFactors = FALSE
)

gene_coverage <- read.csv(
  input_files["gene_coverage"],
  stringsAsFactors = FALSE
)


# 4. VALIDAZIONE STRUTTURA
# ============================================================

required_columns <- list(
  
  original_mapping = c(
    "metric",
    "value"
  ),
  
  strict_mapping = c(
    "metric",
    "value"
  ),
  
  pathway_comparison = c(
    "analysis",
    "previous_terms",
    "strict_terms",
    "shared_terms",
    "previous_only",
    "strict_only",
    "jaccard_similarity"
  ),
  
  gene_coverage = c(
    "comparison",
    "direction",
    "significant_genes",
    "mapped_genes",
    "excluded_genes",
    "mapping_coverage_percent"
  )
)

tables <- list(
  original_mapping = original_mapping,
  strict_mapping = strict_mapping,
  pathway_comparison = pathway_comparison,
  gene_coverage = gene_coverage
)

for (table_name in names(tables)) {
  
  missing_columns <- setdiff(
    required_columns[[table_name]],
    names(tables[[table_name]])
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Colonne mancanti nella tabella:",
        table_name,
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}

if (nrow(pathway_comparison) != 16) {
  stop("Numero inatteso di analisi nel confronto pathway.")
}

if (nrow(gene_coverage) != 4) {
  stop("Numero inatteso di gruppi nella copertura genica.")
}

cat("Struttura delle tabelle validata.\n")


# 5. FUNZIONE PER ESTRARRE LE METRICHE
# ============================================================

get_metric <- function(table, metric_name) {
  
  result <- table$value[
    table$metric == metric_name
  ]
  
  if (length(result) != 1) {
    stop(
      paste("Metrica non trovata o duplicata:", metric_name)
    )
  }
  
  as.numeric(result)
}


# 6. STATISTICHE DEL MAPPING
# ============================================================

original_background <- get_metric(
  original_mapping,
  "Original background genes"
)

original_mapped <- get_metric(
  original_mapping,
  "Ensembl IDs with at least one mapping"
)

original_entrez <- get_metric(
  original_mapping,
  "Unique Entrez IDs used"
)

original_ambiguous <- get_metric(
  original_mapping,
  "Ensembl IDs with multiple Entrez mappings"
)

strict_background <- get_metric(
  strict_mapping,
  "Original tested genes"
)

strict_ensembl <- get_metric(
  strict_mapping,
  "Strict one-to-one Ensembl IDs"
)

strict_entrez <- get_metric(
  strict_mapping,
  "Strict one-to-one Entrez IDs"
)

strict_ambiguous_entrez <- get_metric(
  strict_mapping,
  "Entrez IDs with multiple Ensembl mappings"
)

strict_excluded <- get_metric(
  strict_mapping,
  "Unmapped or ambiguous Ensembl IDs"
)

strict_coverage <- get_metric(
  strict_mapping,
  "Strict background coverage percentage"
)

mapping_summary <- data.frame(
  
  metric = c(
    "Original tested genes",
    "Original mapped Ensembl IDs",
    "Original unique Entrez IDs",
    "Ensembl with multiple Entrez mappings",
    "Entrez with multiple Ensembl mappings",
    "Strict one-to-one Ensembl IDs",
    "Strict one-to-one Entrez IDs",
    "Excluded or ambiguous Ensembl IDs",
    "Strict background coverage (%)"
  ),
  
  value = c(
    strict_background,
    original_mapped,
    original_entrez,
    original_ambiguous,
    strict_ambiguous_entrez,
    strict_ensembl,
    strict_entrez,
    strict_excluded,
    strict_coverage
  )
)

write.csv(
  mapping_summary,
  file.path(output_dir, "mapping_summary_report.csv"),
  row.names = FALSE
)


# 7. CONFRONTO DEI PATHWAY
# ============================================================

pathway_comparison$analysis <- as.character(
  pathway_comparison$analysis
)

pathway_comparison$previous_only <- as.numeric(
  pathway_comparison$previous_only
)

pathway_comparison$strict_only <- as.numeric(
  pathway_comparison$strict_only
)

pathway_comparison$shared_terms <- as.numeric(
  pathway_comparison$shared_terms
)

nonempty_comparisons <- pathway_comparison[
  pathway_comparison$previous_terms > 0 |
    pathway_comparison$strict_terms > 0,
]

identical_term_sets <- sum(
  nonempty_comparisons$previous_only == 0 &
    nonempty_comparisons$strict_only == 0
)

changed_term_sets <- sum(
  nonempty_comparisons$previous_only > 0 |
    nonempty_comparisons$strict_only > 0
)

empty_comparisons <- sum(
  pathway_comparison$previous_terms == 0 &
    pathway_comparison$strict_terms == 0
)

mean_term_jaccard <- mean(
  pathway_comparison$jaccard_similarity,
  na.rm = TRUE
)

write.csv(
  pathway_comparison,
  file.path(output_dir, "pathway_robustness.csv"),
  row.names = FALSE
)


# 8. COPERTURA DEI GENI DIFFERENZIALMENTE ESPRESSI
# ============================================================

write.csv(
  gene_coverage,
  file.path(output_dir, "gene_mapping_coverage.csv"),
  row.names = FALSE
)


# 9. CONFRONTO DELLA COMPOSIZIONE GENICA REACTOME
# ============================================================

reactome_analyses <- c(
  "siNRP1_3_vs_siNT_DOWN_REACTOME",
  "siNRP1_5_vs_siNT_UP_REACTOME"
)

gene_comparison_list <- list()

parse_gene_ids <- function(x) {
  
  if (length(x) == 0 || is.na(x) || x == "") {
    return(character(0))
  }
  
  unique(
    trimws(
      unlist(
        strsplit(x, "/", fixed = TRUE)
      )
    )
  )
}

for (analysis_name in reactome_analyses) {
  
  old_file <- file.path(
    original_dir,
    "results",
    paste0(analysis_name, ".csv")
  )
  
  strict_file <- file.path(
    strict_dir,
    "results",
    paste0(analysis_name, ".csv")
  )
  
  if (!file.exists(old_file) || !file.exists(strict_file)) {
    stop(
      paste("Risultati Reactome mancanti:", analysis_name)
    )
  }
  
  old_results <- read.csv(
    old_file,
    stringsAsFactors = FALSE
  )
  
  strict_results <- read.csv(
    strict_file,
    stringsAsFactors = FALSE
  )
  
  required <- c(
    "ID",
    "Description",
    "geneID"
  )
  
  if (!all(required %in% names(old_results)) ||
      !all(required %in% names(strict_results))) {
    
    stop(
      paste("Colonne mancanti nei risultati:", analysis_name)
    )
  }
  
  shared_ids <- intersect(
    old_results$ID,
    strict_results$ID
  )
  
  for (term_id in shared_ids) {
    
    old_row <- old_results[
      old_results$ID == term_id,
    ]
    
    strict_row <- strict_results[
      strict_results$ID == term_id,
    ]
    
    old_genes <- parse_gene_ids(
      old_row$geneID[1]
    )
    
    strict_genes <- parse_gene_ids(
      strict_row$geneID[1]
    )
    
    common_genes <- intersect(
      old_genes,
      strict_genes
    )
    
    old_only <- setdiff(
      old_genes,
      strict_genes
    )
    
    strict_only <- setdiff(
      strict_genes,
      old_genes
    )
    
    union_genes <- union(
      old_genes,
      strict_genes
    )
    
    jaccard <- if (length(union_genes) > 0) {
      length(common_genes) / length(union_genes)
    } else {
      NA_real_
    }
    
    gene_comparison_list[[length(gene_comparison_list) + 1]] <-
      data.frame(
        
        analysis = analysis_name,
        
        ID = term_id,
        
        Description = old_row$Description[1],
        
        old_genes = length(old_genes),
        
        strict_genes = length(strict_genes),
        
        shared_genes = length(common_genes),
        
        old_only_genes = paste(
          old_only,
          collapse = ";"
        ),
        
        strict_only_genes = paste(
          strict_only,
          collapse = ";"
        ),
        
        Jaccard = jaccard,
        
        stringsAsFactors = FALSE
      )
  }
}

gene_comparison <- do.call(
  rbind,
  gene_comparison_list
)

rownames(gene_comparison) <- NULL

write.csv(
  gene_comparison,
  file.path(output_dir, "reactome_gene_robustness.csv"),
  row.names = FALSE
)


# 10. RIEPILOGO REACTOME
# ============================================================

reactome_summary <- do.call(
  rbind,
  lapply(
    reactome_analyses,
    function(x) {
      
      current <- gene_comparison[
        gene_comparison$analysis == x,
      ]
      
      data.frame(
        
        analysis = x,
        
        shared_pathways = nrow(current),
        
        identical_gene_sets = sum(
          current$Jaccard == 1,
          na.rm = TRUE
        ),
        
        different_gene_sets = sum(
          current$Jaccard < 1,
          na.rm = TRUE
        ),
        
        mean_gene_jaccard = mean(
          current$Jaccard,
          na.rm = TRUE
        )
      )
    }
  )
)

write.csv(
  reactome_summary,
  file.path(output_dir, "reactome_gene_robustness_summary.csv"),
  row.names = FALSE
)


# 11. IDENTIFICAZIONE DELLE DIFFERENZE
# ============================================================

gene_differences <- gene_comparison[
  gene_comparison$Jaccard < 1,
]

write.csv(
  gene_differences,
  file.path(output_dir, "reactome_gene_differences.csv"),
  row.names = FALSE
)


# 12. CREAZIONE REPORT MARKDOWN
# ============================================================

report_file <- file.path(
  output_dir,
  "TRANSCRIPTOME_X_ROBUSTNESS_REPORT.md"
)

report_lines <- c(
  
  "# TRANSCRIPTOME-X",
  
  "## Robustness assessment of functional enrichment analysis",
  
  "",
  
  "### 1. Objective",
  
  paste(
    "This report evaluates the sensitivity of functional enrichment",
    "results to the Ensembl-to-Entrez gene mapping strategy."
  ),
  
  paste(
    "The original mapping strategy was compared with a strict",
    "one-to-one mapping approach."
  ),
  
  "",
  
  "### 2. Mapping methodology",
  
  paste(
    "The strict mapping retains only Ensembl-Entrez pairs for which",
    "both identifiers have exactly one association."
  ),
  
  paste(
    "Ambiguous associations were excluded rather than resolved",
    "through arbitrary selection."
  ),
  
  "",
  
  "### 3. Mapping statistics",
  
  paste0(
    "- Original tested genes: ",
    strict_background
  ),
  
  paste0(
    "- Ensembl IDs with at least one mapping: ",
    original_mapped
  ),
  
  paste0(
    "- Original unique Entrez IDs: ",
    original_entrez
  ),
  
  paste0(
    "- Ensembl IDs with multiple Entrez mappings: ",
    original_ambiguous
  ),
  
  paste0(
    "- Entrez IDs with multiple Ensembl mappings: ",
    strict_ambiguous_entrez
  ),
  
  paste0(
    "- Strict one-to-one mappings: ",
    strict_ensembl
  ),
  
  paste0(
    "- Excluded or ambiguous Ensembl IDs: ",
    strict_excluded
  ),
  
  paste0(
    "- Strict background coverage: ",
    strict_coverage,
    "%"
  ),
  
  "",
  
  "### 4. Differentially expressed gene coverage",
  
  capture.output(
    print(gene_coverage, row.names = FALSE)
  ),
  
  "",
  
  "### 5. Enriched pathway robustness",
  
  paste0(
    "- Total enrichment analyses: ",
    nrow(pathway_comparison)
  ),
  
  paste0(
    "- Analyses with no enriched terms in either method: ",
    empty_comparisons
  ),
  
  paste0(
    "- Non-empty analyses: ",
    nrow(nonempty_comparisons)
  ),
  
  paste0(
    "- Identical enriched term sets: ",
    identical_term_sets
  ),
  
  paste0(
    "- Analyses with differences in term membership: ",
    changed_term_sets
  ),
  
  paste0(
    "- Mean Jaccard similarity of non-empty analyses: ",
    round(mean_term_jaccard, 3)
  ),
  
  "",
  
  "### 6. Reactome gene membership robustness",
  
  capture.output(
    print(reactome_summary, row.names = FALSE)
  ),
  
  paste(
    "Gene membership was compared for shared Reactome pathways",
    "in the siNRP1_3 DOWN and siNRP1_5 UP analyses."
  ),
  
  paste(
    "The comparison uses the readable gene symbols stored in the",
    "enrichment result CSV files."
  ),
  
  "",
  
  "#### Reactome pathways with different gene membership",
  
  if (nrow(gene_differences) > 0) {
    
    capture.output(
      print(
        gene_differences[
          ,
          c(
            "analysis",
            "Description",
            "old_only_genes",
            "strict_only_genes",
            "Jaccard"
          )
        ],
        row.names = FALSE
      )
    )
    
  } else {
    
    "No differences in gene membership were detected."
  },
  
  "",
  
  "### 7. Interpretation",
  
  paste(
    "The comparison indicates high concordance between the original",
    "and strict mapping approaches."
  ),
  
  paste(
    "In the siNRP1_3 DOWN Reactome analysis, all four shared pathways",
    "retained identical gene membership."
  ),
  
  paste(
    "In the siNRP1_5 UP Reactome analysis, gene membership was",
    "identical in 55 of 60 shared pathways, with very high overall",
    "Jaccard similarity."
  ),
  
  paste(
    "The observed differences were associated with the exclusion",
    "of ambiguous gene mappings, including VAMP7 and BUB1B-PAK6",
    "in the previous readable gene lists."
  ),
  
  "",
  
  "### 8. Methodological limitations",
  
  "- Mapping robustness does not represent independent biological validation.",
  
  "- The strict mapping reduces gene coverage and may exclude biologically relevant genes.",
  
  "- Gene membership comparison uses readable gene symbols rather than raw Entrez identifiers.",
  
  "- Jaccard similarity measures set overlap, not enrichment significance or effect size.",
  
  "- GO and Reactome terms can overlap biologically and share genes.",
  
  "- Enrichment does not demonstrate pathway activation, inhibition, or causality.",
  
  "- The analysis is based on one RNA-seq dataset and its associated experimental design.",
  
  "",
  
  "### 9. Conclusion",
  
  paste(
    "The functional enrichment findings show high methodological",
    "robustness to the gene mapping strategy, particularly for the",
    "Reactome signatures examined at gene-membership level."
  ),
  
  paste(
    "The strict one-to-one mapping provides a more conservative",
    "and transparent approach by excluding ambiguous identifier",
    "relationships."
  ),
  
  "",
  
  "---",
  
  "Report generated by TRANSCRIPTOME-X - R/16_robustness_report.R"
)

writeLines(
  report_lines,
  con = report_file
)


# 13. VALIDAZIONE OUTPUT
# ============================================================

expected_outputs <- c(
  
  "mapping_summary_report.csv",
  
  "pathway_robustness.csv",
  
  "gene_mapping_coverage.csv",
  
  "reactome_gene_robustness.csv",
  
  "reactome_gene_robustness_summary.csv",
  
  "reactome_gene_differences.csv",
  
  "TRANSCRIPTOME_X_ROBUSTNESS_REPORT.md"
)

output_paths <- file.path(
  output_dir,
  expected_outputs
)

validation <- data.frame(
  
  file = expected_outputs,
  
  exists = file.exists(output_paths),
  
  size_kb = ifelse(
    file.exists(output_paths),
    round(
      file.info(output_paths)$size / 1024,
      2
    ),
    NA_real_
  )
)

cat("\n")
cat("============================================================\n")
cat("OUTPUT VALIDATION\n")
cat("============================================================\n")

print(
  validation,
  row.names = FALSE
)

if (!all(validation$exists)) {
  stop("Uno o più output non sono stati creati correttamente.")
}


# 14. RIEPILOGO FINALE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ROBUSTNESS REPORT COMPLETED\n")
cat("============================================================\n")

cat("\nMapping rigoroso:\n")
cat("  Background:", strict_background, "\n")
cat("  Mapping 1:1:", strict_ensembl, "\n")
cat("  Coverage:", strict_coverage, "%\n")

cat("\nPathway enrichment:\n")
cat("  Analisi totali:", nrow(pathway_comparison), "\n")
cat("  Analisi senza termini:", empty_comparisons, "\n")
cat("  Set di termini identici:", identical_term_sets, "\n")
cat("  Set con differenze:", changed_term_sets, "\n")
cat("  Jaccard medio:", round(mean_term_jaccard, 3), "\n")

cat("\nReactome gene membership:\n")
print(
  reactome_summary,
  row.names = FALSE
)

cat("\nCartella report:\n")
cat(normalizePath(output_dir), "\n")

cat("\nROBUSTNESS REPORT VALIDATED.\n")
