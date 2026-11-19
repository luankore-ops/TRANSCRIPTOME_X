
# ============================================================
# TRANSCRIPTOME-X
# 20 - CORE GENE ANALYSIS (ENSEMBL-BASED)
# ============================================================
#
# Dataset: GSE266566
# Cell model: MDA-MB-231
#
# Changes from script 18:
# - Ensembl stable ID is the integration key.
# - Enrichment and DE gene symbols are kept separately.
# - Ambiguous symbol mappings are not resolved arbitrarily.
# - Original outputs are preserved in the old directory.
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - CORE GENE ANALYSIS V2\n")
cat("ENSEMBL-BASED INTEGRATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

project_root <- getwd()

input_dir <- file.path(
  project_root,
  "data/processed/functional_enrichment_strict/reports",
  "integrated_biological_interpretation"
)

de_dir <- file.path(
  project_root,
  "data/processed/differential_expression/annotated"
)

output_dir <- file.path(
  project_root,
  "data/processed/functional_enrichment_strict/reports",
  "core_gene_analysis_v2"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

comparisons <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

directions <- c("UP", "DOWN")


# ------------------------------------------------------------
# 2. CHECK PACKAGES
# ------------------------------------------------------------

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("Pacchetto AnnotationDbi non installato.")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Pacchetto org.Hs.eg.db non installato.")
}

annotation_db <- org.Hs.eg.db::org.Hs.eg.db

cat(
  "AnnotationDbi:",
  as.character(packageVersion("AnnotationDbi")),
  "\n"
)

cat(
  "org.Hs.eg.db:",
  as.character(packageVersion("org.Hs.eg.db")),
  "\n"
)


# ------------------------------------------------------------
# 3. LOAD ENRICHMENT RESULTS
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE ENRICHMENT RESULTS\n")
cat("------------------------------------------------------------\n")

enrichment_file <- file.path(
  input_dir,
  "integrated_enrichment_results.csv"
)

if (!file.exists(enrichment_file)) {
  stop("File enrichment non trovato: ", enrichment_file)
}

enrichment <- read.csv(
  enrichment_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_enrichment_columns <- c(
  "ID",
  "Description",
  "p.adjust",
  "geneID",
  "comparison",
  "direction",
  "method"
)

missing_columns <- setdiff(
  required_enrichment_columns,
  colnames(enrichment)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonne mancanti nel file enrichment: ",
    paste(missing_columns, collapse = ", ")
  )
}

enrichment <- enrichment[
  !is.na(enrichment$geneID) &
    enrichment$geneID != "",
  ,
  drop = FALSE
]

cat("Termini con geni associati:", nrow(enrichment), "\n")
cat(
  "Simboli distinti:",
  length(unique(unlist(strsplit(
    enrichment$geneID,
    "/",
    fixed = TRUE
  )))),
  "\n"
)


# ------------------------------------------------------------
# 4. LOAD DE RESULTS
# ------------------------------------------------------------

cat("\n2. IMPORTAZIONE DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")

de_list <- list()

for (comparison in comparisons) {
  
  de_file <- file.path(
    de_dir,
    paste0(comparison, "_annotated.csv")
  )
  
  if (!file.exists(de_file)) {
    stop("File DE non trovato: ", de_file)
  }
  
  de <- read.csv(
    de_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_de_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  missing_de_columns <- setdiff(
    required_de_columns,
    colnames(de)
  )
  
  if (length(missing_de_columns) > 0) {
    stop(
      "Colonne mancanti in ",
      de_file,
      ": ",
      paste(missing_de_columns, collapse = ", ")
    )
  }
  
  # Preserve the original DE symbol.
  de$de_gene_name <- de$gene_name
  
  # Remove Ensembl version suffix:
  # ENSG00000188522.16 -> ENSG00000188522
  de$ensembl_id <- sub(
    "\\..*$",
    "",
    as.character(de$gene_id)
  )
  
  de$comparison <- comparison
  
  de_list[[comparison]] <- de
}

de_all <- do.call(rbind, de_list)

rownames(de_all) <- NULL

de_all <- de_all[
  !is.na(de_all$ensembl_id) &
    de_all$ensembl_id != "",
  ,
  drop = FALSE
]

# A stable Ensembl ID must identify one DE row per comparison.
de_duplicate_key <- duplicated(
  de_all[, c("comparison", "ensembl_id")]
)

if (any(de_duplicate_key)) {
  duplicated_rows <- de_all[
    de_duplicate_key,
    c("comparison", "gene_id", "gene_name", "ensembl_id"),
    drop = FALSE
  ]
  
  write.csv(
    duplicated_rows,
    file.path(output_dir, "DE_duplicate_ensembl_ids.csv"),
    row.names = FALSE
  )
  
  stop(
    "Trovati identificativi Ensembl duplicati nei file DE. ",
    "Controllare DE_duplicate_ensembl_ids.csv."
  )
}

cat("Righe DE caricate:", nrow(de_all), "\n")
cat(
  "Identificativi Ensembl distinti:",
  length(unique(de_all$ensembl_id)),
  "\n"
)


# ------------------------------------------------------------
# 5. BUILD SYMBOL -> ENSEMBL MAPPING
# ------------------------------------------------------------

cat("\n3. MAPPING SIMBOLO-ENSEMBL\n")
cat("------------------------------------------------------------\n")

enrichment_symbols <- unique(
  trimws(
    unlist(
      strsplit(
        enrichment$geneID,
        "/",
        fixed = TRUE
      )
    )
  )
)

enrichment_symbols <- sort(
  enrichment_symbols[
    !is.na(enrichment_symbols) &
      enrichment_symbols != ""
  ]
)

symbol_mapping_raw <- AnnotationDbi::select(
  annotation_db,
  keys = enrichment_symbols,
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENSEMBL",
    "ENTREZID"
  )
)

symbol_mapping_raw <- unique(
  symbol_mapping_raw[
    !is.na(symbol_mapping_raw$SYMBOL) &
      !is.na(symbol_mapping_raw$ENSEMBL) &
      !is.na(symbol_mapping_raw$ENTREZID),
    ,
    drop = FALSE
  ]
)

# Count distinct Ensembl IDs per symbol.
ensembl_per_symbol <- split(
  symbol_mapping_raw$ENSEMBL,
  symbol_mapping_raw$SYMBOL
)

symbol_ensembl_counts <- data.frame(
  SYMBOL = names(ensembl_per_symbol),
  n_ensembl = vapply(
    ensembl_per_symbol,
    function(x) length(unique(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

# Count distinct Entrez IDs per symbol.
entrez_per_symbol <- split(
  symbol_mapping_raw$ENTREZID,
  symbol_mapping_raw$SYMBOL
)

symbol_entrez_counts <- data.frame(
  SYMBOL = names(entrez_per_symbol),
  n_entrez = vapply(
    entrez_per_symbol,
    function(x) length(unique(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

# Count symbols per Ensembl ID.
symbols_per_ensembl <- split(
  symbol_mapping_raw$SYMBOL,
  symbol_mapping_raw$ENSEMBL
)

ensembl_symbol_counts <- data.frame(
  ENSEMBL = names(symbols_per_ensembl),
  n_symbols = vapply(
    symbols_per_ensembl,
    function(x) length(unique(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

# Count symbols per Entrez ID.
symbols_per_entrez <- split(
  symbol_mapping_raw$SYMBOL,
  symbol_mapping_raw$ENTREZID
)

entrez_symbol_counts <- data.frame(
  ENTREZID = names(symbols_per_entrez),
  n_symbols = vapply(
    symbols_per_entrez,
    function(x) length(unique(x)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

# Combine mapping diagnostics.
symbol_audit <- data.frame(
  SYMBOL = enrichment_symbols,
  stringsAsFactors = FALSE
)

symbol_audit <- merge(
  symbol_audit,
  symbol_ensembl_counts,
  by = "SYMBOL",
  all.x = TRUE,
  sort = FALSE
)

symbol_audit <- merge(
  symbol_audit,
  symbol_entrez_counts,
  by = "SYMBOL",
  all.x = TRUE,
  sort = FALSE
)

symbol_audit$n_ensembl[
  is.na(symbol_audit$n_ensembl)
] <- 0L

symbol_audit$n_entrez[
  is.na(symbol_audit$n_entrez)
] <- 0L

symbol_audit$mapping_status <- ifelse(
  symbol_audit$n_ensembl == 0 |
    symbol_audit$n_entrez == 0,
  "NO_MAPPING",
  ifelse(
    symbol_audit$n_ensembl != 1 |
      symbol_audit$n_entrez != 1,
    "AMBIGUOUS_SYMBOL",
    "CANDIDATE_UNIQUE"
  )
)

# Keep only symbols with one Ensembl and one Entrez.
candidate_mapping <- symbol_mapping_raw[
  symbol_audit$SYMBOL[
    match(
      symbol_mapping_raw$SYMBOL,
      symbol_audit$SYMBOL
    )
  ] %in% symbol_audit$SYMBOL[
    symbol_audit$mapping_status == "CANDIDATE_UNIQUE"
  ],
  ,
  drop = FALSE
]

# Also require each Ensembl and Entrez ID to map to one symbol.
candidate_mapping <- candidate_mapping[
  candidate_mapping$ENSEMBL %in%
    ensembl_symbol_counts$ENSEMBL[
      ensembl_symbol_counts$n_symbols == 1
    ] &
    candidate_mapping$ENTREZID %in%
    entrez_symbol_counts$ENTREZID[
      entrez_symbol_counts$n_symbols == 1
    ],
  ,
  drop = FALSE
]

strict_symbol_mapping <- unique(
  candidate_mapping[
    ,
    c("SYMBOL", "ENSEMBL", "ENTREZID"),
    drop = FALSE
  ]
)

# Final defensive check.
if (anyDuplicated(strict_symbol_mapping$SYMBOL) ||
    anyDuplicated(strict_symbol_mapping$ENSEMBL) ||
    anyDuplicated(strict_symbol_mapping$ENTREZID)) {
  stop("Il mapping strict non è univoco.")
}

symbol_audit$mapping_status[
  symbol_audit$SYMBOL %in% strict_symbol_mapping$SYMBOL
] <- "STRICT_ONE_TO_ONE"

symbol_audit$mapping_status[
  symbol_audit$mapping_status == "CANDIDATE_UNIQUE"
] <- "AMBIGUOUS_REVERSE_MAPPING"

write.csv(
  symbol_audit,
  file.path(output_dir, "gene_symbol_mapping_audit.csv"),
  row.names = FALSE
)

write.csv(
  strict_symbol_mapping,
  file.path(output_dir, "gene_symbol_mapping_strict.csv"),
  row.names = FALSE
)

cat("Simboli enrichment:", length(enrichment_symbols), "\n")
cat(
  "Mapping strict 1:1:",
  nrow(strict_symbol_mapping),
  "\n"
)
cat(
  "Simboli senza mapping:",
  sum(symbol_audit$mapping_status == "NO_MAPPING"),
  "\n"
)
cat(
  "Simboli non univoci o ambigui:",
  sum(
    symbol_audit$mapping_status %in% c(
      "AMBIGUOUS_SYMBOL",
      "AMBIGUOUS_REVERSE_MAPPING"
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. EXPAND GENE MEMBERSHIP WITH STABLE IDS
# ------------------------------------------------------------

cat("\n4. ESPANSIONE GENE MEMBERSHIP\n")
cat("------------------------------------------------------------\n")

gene_rows <- vector("list", nrow(enrichment))

for (i in seq_len(nrow(enrichment))) {
  
  genes <- unique(
    trimws(
      unlist(
        strsplit(
          enrichment$geneID[i],
          "/",
          fixed = TRUE
        )
      )
    )
  )
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  if (length(genes) == 0) {
    next
  }
  
  map_index <- match(
    genes,
    strict_symbol_mapping$SYMBOL
  )
  
  status_index <- match(
    genes,
    symbol_audit$SYMBOL
  )
  
  current <- data.frame(
    comparison = enrichment$comparison[i],
    direction = enrichment$direction[i],
    method = enrichment$method[i],
    term_id = enrichment$ID[i],
    term_description = enrichment$Description[i],
    adjusted_p_value = enrichment$`p.adjust`[i],
    enrichment_symbol = genes,
    ensembl_id = strict_symbol_mapping$ENSEMBL[map_index],
    entrez_id = strict_symbol_mapping$ENTREZID[map_index],
    mapping_status = symbol_audit$mapping_status[status_index],
    stringsAsFactors = FALSE
  )
  
  # Compatibility field: this is the enrichment symbol,
  # not the DE annotation symbol.
  current$gene_name <- current$enrichment_symbol
  
  gene_rows[[i]] <- current
}

gene_rows <- gene_rows[
  !vapply(gene_rows, is.null, logical(1))
]

if (length(gene_rows) == 0) {
  stop("Nessun gene estratto dai risultati enrichment.")
}

gene_membership <- do.call(
  rbind,
  gene_rows
)

rownames(gene_membership) <- NULL

gene_membership <- unique(gene_membership)

# A stable key is used for recurrence:
# Ensembl ID when available; otherwise a symbol-specific key.
gene_membership$gene_key <- ifelse(
  !is.na(gene_membership$ensembl_id),
  gene_membership$ensembl_id,
  paste0(
    "UNMAPPED_SYMBOL:",
    gene_membership$enrichment_symbol
  )
)

write.csv(
  gene_membership,
  file.path(output_dir, "gene_pathway_membership_ensembl.csv"),
  row.names = FALSE
)

cat("Associazioni gene-term:", nrow(gene_membership), "\n")
cat(
  "Simboli distinti:",
  length(unique(gene_membership$enrichment_symbol)),
  "\n"
)
cat(
  "Geni con Ensembl:",
  length(unique(
    gene_membership$ensembl_id[
      !is.na(gene_membership$ensembl_id)
    ]
  )),
  "\n"
)
cat(
  "Geni senza mapping strict:",
  length(unique(
    gene_membership$gene_key[
      is.na(gene_membership$ensembl_id)
    ]
  )),
  "\n"
)


# ------------------------------------------------------------
# 7. FREQUENCY BY COMPARISON, DIRECTION AND METHOD
# ------------------------------------------------------------

cat("\n5. FREQUENZA CORE GENES PER METODO\n")
cat("------------------------------------------------------------\n")

method_key <- interaction(
  gene_membership$comparison,
  gene_membership$direction,
  gene_membership$method,
  drop = TRUE,
  lex.order = TRUE
)

method_groups <- split(
  gene_membership,
  method_key
)

frequency_results <- lapply(
  method_groups,
  function(current) {
    
    gene_groups <- split(
      current,
      current$gene_key
    )
    
    do.call(
      rbind,
      lapply(
        gene_groups,
        function(g) {
          
          data.frame(
            comparison = g$comparison[1],
            direction = g$direction[1],
            method = g$method[1],
            gene_key = g$gene_key[1],
            enrichment_symbol = g$enrichment_symbol[1],
            ensembl_id = g$ensembl_id[1],
            entrez_id = g$entrez_id[1],
            mapping_status = g$mapping_status[1],
            pathway_count = length(unique(g$term_id)),
            term_ids = paste(
              unique(g$term_id),
              collapse = ";"
            ),
            term_descriptions = paste(
              unique(g$term_description),
              collapse = " | "
            ),
            best_adjusted_p_value = min(
              g$adjusted_p_value,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }
)

gene_frequency <- do.call(
  rbind,
  frequency_results
)

rownames(gene_frequency) <- NULL

gene_frequency <- gene_frequency[
  order(
    gene_frequency$comparison,
    gene_frequency$direction,
    gene_frequency$method,
    -gene_frequency$pathway_count,
    gene_frequency$best_adjusted_p_value
  ),
  ,
  drop = FALSE
]

write.csv(
  gene_frequency,
  file.path(output_dir, "core_gene_frequency_by_method.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. OVERALL FREQUENCY
# ------------------------------------------------------------

cat("\n6. FREQUENZA COMPLESSIVA\n")
cat("------------------------------------------------------------\n")

overall_key <- interaction(
  gene_membership$comparison,
  gene_membership$direction,
  drop = TRUE,
  lex.order = TRUE
)

overall_groups <- split(
  gene_membership,
  overall_key
)

overall_results <- lapply(
  overall_groups,
  function(current) {
    
    gene_groups <- split(
      current,
      current$gene_key
    )
    
    do.call(
      rbind,
      lapply(
        gene_groups,
        function(g) {
          
          unique_terms <- unique(
            g[, c("method", "term_id")]
          )
          
          data.frame(
            comparison = g$comparison[1],
            direction = g$direction[1],
            gene_key = g$gene_key[1],
            enrichment_symbol = g$enrichment_symbol[1],
            ensembl_id = g$ensembl_id[1],
            entrez_id = g$entrez_id[1],
            mapping_status = g$mapping_status[1],
            total_term_associations = nrow(g),
            unique_terms = nrow(unique_terms),
            methods = paste(
              sort(unique(g$method)),
              collapse = ";"
            ),
            term_ids = paste(
              paste(
                unique_terms$method,
                unique_terms$term_id,
                sep = ":"
              ),
              collapse = ";"
            ),
            best_adjusted_p_value = min(
              g$adjusted_p_value,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }
)

overall_frequency <- do.call(
  rbind,
  overall_results
)

rownames(overall_frequency) <- NULL

overall_frequency <- overall_frequency[
  order(
    overall_frequency$comparison,
    overall_frequency$direction,
    -overall_frequency$unique_terms,
    overall_frequency$best_adjusted_p_value
  ),
  ,
  drop = FALSE
]

write.csv(
  overall_frequency,
  file.path(output_dir, "core_gene_frequency_overall.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. INTEGRATE BY ENSEMBL ID
# ------------------------------------------------------------

cat("\n7. INTEGRAZIONE DE CON IDENTIFICATIVO ENSEMBL\n")
cat("------------------------------------------------------------\n")

de_mapping <- unique(
  de_all[
    ,
    c(
      "comparison",
      "ensembl_id",
      "gene_id",
      "de_gene_name",
      "logFC",
      "FDR"
    ),
    drop = FALSE
  ]
)

# Only mapped genes are merged by Ensembl.
mapped_core <- overall_frequency[
  !is.na(overall_frequency$ensembl_id),
  ,
  drop = FALSE
]

unmapped_core <- overall_frequency[
  is.na(overall_frequency$ensembl_id),
  ,
  drop = FALSE
]

integrated_mapped <- merge(
  mapped_core,
  de_mapping,
  by = c("comparison", "ensembl_id"),
  all.x = TRUE,
  sort = FALSE
)

# Preserve unmapped enrichment genes without matching
# unrelated DE rows that also have missing IDs.
unmapped_core$gene_id <- NA_character_
unmapped_core$de_gene_name <- NA_character_
unmapped_core$logFC <- NA_real_
unmapped_core$FDR <- NA_real_

integrated_core_genes <- rbind(
  integrated_mapped,
  unmapped_core
)

rownames(integrated_core_genes) <- NULL

integrated_core_genes$de_significant <- (
  !is.na(integrated_core_genes$FDR) &
    integrated_core_genes$FDR < 0.05
)

integrated_core_genes$de_direction <- ifelse(
  is.na(integrated_core_genes$logFC),
  NA_character_,
  ifelse(
    integrated_core_genes$logFC > 0,
    "UP",
    ifelse(
      integrated_core_genes$logFC < 0,
      "DOWN",
      "ZERO"
    )
  )
)

integrated_core_genes$direction_consistent <- (
  integrated_core_genes$de_significant &
    integrated_core_genes$direction ==
    integrated_core_genes$de_direction
)

integrated_core_genes$de_mapping_status <- ifelse(
  is.na(integrated_core_genes$ensembl_id),
  "NO_UNAMBIGUOUS_ENSEMBL_MAPPING",
  ifelse(
    is.na(integrated_core_genes$gene_id),
    "ENSEMBL_NOT_FOUND_IN_DE",
    "MATCHED_BY_ENSEMBL"
  )
)

integrated_core_genes <- integrated_core_genes[
  order(
    integrated_core_genes$comparison,
    integrated_core_genes$direction,
    -integrated_core_genes$unique_terms,
    integrated_core_genes$FDR,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  integrated_core_genes,
  file.path(
    output_dir,
    "core_genes_with_differential_expression.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. TOP 20 CORE GENES
# ------------------------------------------------------------

cat("\n8. TOP CORE GENES\n")
cat("------------------------------------------------------------\n")

top_core_list <- list()

for (comparison in comparisons) {
  
  for (direction in directions) {
    
    current <- integrated_core_genes[
      integrated_core_genes$comparison == comparison &
        integrated_core_genes$direction == direction,
      ,
      drop = FALSE
    ]
    
    current <- current[
      order(
        -current$unique_terms,
        current$FDR,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    if (nrow(current) > 0) {
      top_core_list[[
        paste(comparison, direction, sep = "_")
      ]] <- head(current, 20)
    }
  }
}

if (length(top_core_list) > 0) {
  
  top_core_genes <- do.call(
    rbind,
    top_core_list
  )
  
  rownames(top_core_genes) <- NULL
  
} else {
  
  top_core_genes <- integrated_core_genes[
    0,
    ,
    drop = FALSE
  ]
}

write.csv(
  top_core_genes,
  file.path(output_dir, "top_20_core_genes.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. MAPPING AND INTEGRATION AUDIT
# ------------------------------------------------------------

cat("\n9. AUDIT MAPPING E INTEGRAZIONE\n")
cat("------------------------------------------------------------\n")

mapping_summary <- as.data.frame(
  table(symbol_audit$mapping_status),
  stringsAsFactors = FALSE
)

names(mapping_summary) <- c(
  "mapping_status",
  "symbol_count"
)

write.csv(
  mapping_summary,
  file.path(output_dir, "mapping_status_summary.csv"),
  row.names = FALSE
)

integration_audit <- as.data.frame(
  table(integrated_core_genes$de_mapping_status),
  stringsAsFactors = FALSE
)

names(integration_audit) <- c(
  "de_mapping_status",
  "gene_count"
)

write.csv(
  integration_audit,
  file.path(output_dir, "DE_integration_status_summary.csv"),
  row.names = FALSE
)

cat("\nMapping status:\n")
print(mapping_summary, row.names = FALSE)

cat("\nDE integration status:\n")
print(integration_audit, row.names = FALSE)


# ------------------------------------------------------------
# 12. SUMMARY
# ------------------------------------------------------------

cat("\n10. RIEPILOGO CORE GENE ANALYSIS\n")
cat("------------------------------------------------------------\n")

summary_table <- do.call(
  rbind,
  lapply(
    comparisons,
    function(comparison) {
      
      do.call(
        rbind,
        lapply(
          directions,
          function(direction) {
            
            current <- integrated_core_genes[
              integrated_core_genes$comparison == comparison &
                integrated_core_genes$direction == direction,
              ,
              drop = FALSE
            ]
            
            data.frame(
              comparison = comparison,
              direction = direction,
              core_genes = length(unique(current$gene_key)),
              genes_with_ensembl = length(unique(
                current$ensembl_id[
                  !is.na(current$ensembl_id)
                ]
              )),
              genes_matched_to_DE = length(unique(
                current$ensembl_id[
                  !is.na(current$gene_id)
                ]
              )),
              significant_DE_genes = length(unique(
                current$gene_key[
                  current$de_significant
                ]
              )),
              direction_consistent_genes = length(unique(
                current$gene_key[
                  current$direction_consistent
                ]
              )),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

write.csv(
  summary_table,
  file.path(output_dir, "core_gene_summary.csv"),
  row.names = FALSE
)

print(summary_table, row.names = FALSE)


# ------------------------------------------------------------
# 13. MARKDOWN REPORT
# ------------------------------------------------------------

cat("\n11. GENERAZIONE REPORT\n")
cat("------------------------------------------------------------\n")

report <- c(
  "# TRANSCRIPTOME-X",
  "",
  "# Core gene analysis - Ensembl-based integration",
  "",
  "## 1. Methodology",
  "",
  "All significant GO and Reactome enrichment terms were used.",
  "Genes were extracted from the readable geneID field.",
  "Enrichment symbols were mapped to Ensembl stable identifiers.",
  "Only unambiguous one-to-one symbol-Ensembl-Entrez mappings were retained.",
  "Differential expression results were integrated using Ensembl stable IDs.",
  "Enrichment and DE gene symbols were preserved in separate columns.",
  "",
  "## 2. Annotation versions",
  "",
  paste(
    "AnnotationDbi:",
    as.character(packageVersion("AnnotationDbi"))
  ),
  paste(
    "org.Hs.eg.db:",
    as.character(packageVersion("org.Hs.eg.db"))
  ),
  "",
  "## 3. Summary",
  "",
  paste(
    capture.output(
      print(summary_table, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 4. Mapping status",
  "",
  paste(
    capture.output(
      print(mapping_summary, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 5. DE integration status",
  "",
  paste(
    capture.output(
      print(integration_audit, row.names = FALSE)
    ),
    collapse = "\n"
  ),
  "",
  "## 6. Top 20 core genes",
  ""
)

for (comparison in comparisons) {
  
  for (direction in directions) {
    
    current <- top_core_genes[
      top_core_genes$comparison == comparison &
        top_core_genes$direction == direction,
      ,
      drop = FALSE
    ]
    
    report <- c(
      report,
      paste0("### ", comparison, " - ", direction),
      ""
    )
    
    if (nrow(current) == 0) {
      
      report <- c(
        report,
        "No core genes identified.",
        ""
      )
      
    } else {
      
      report_columns <- intersect(
        c(
          "enrichment_symbol",
          "de_gene_name",
          "ensembl_id",
          "unique_terms",
          "methods",
          "logFC",
          "FDR",
          "de_direction",
          "direction_consistent",
          "de_mapping_status"
        ),
        colnames(current)
      )
      
      report_table <- current[
        ,
        report_columns,
        drop = FALSE
      ]
      
      report <- c(
        report,
        paste(
          capture.output(
            print(report_table, row.names = FALSE)
          ),
          collapse = "\n"
        ),
        ""
      )
    }
  }
}

report <- c(
  report,
  "## 7. Important limitations",
  "",
  "- Gene recurrence is not an independent measure of biological importance.",
  "- GO terms are hierarchical and may be redundant.",
  "- Reactome pathways may share many genes.",
  "- Ambiguous symbol mappings are not assigned arbitrarily.",
  "- Genes without an unambiguous mapping remain in recurrence tables but are not integrated with DE.",
  "- The same dataset is used for differential expression and enrichment.",
  "- Recurrence does not establish causality or pathway activation.",
  "- Expression direction is checked against the original DE results.",
  "",
  "## 8. Output files",
  "",
  "- gene_pathway_membership_ensembl.csv",
  "- gene_symbol_mapping_audit.csv",
  "- gene_symbol_mapping_strict.csv",
  "- core_gene_frequency_by_method.csv",
  "- core_gene_frequency_overall.csv",
  "- core_genes_with_differential_expression.csv",
  "- top_20_core_genes.csv",
  "- mapping_status_summary.csv",
  "- DE_integration_status_summary.csv",
  "- core_gene_summary.csv"
)

report_file <- file.path(
  output_dir,
  "TRANSCRIPTOME_X_CORE_GENE_REPORT.md"
)

writeLines(
  report,
  report_file
)


# ------------------------------------------------------------
# 14. OUTPUT VALIDATION
# ------------------------------------------------------------

cat("\n12. OUTPUT VALIDATION\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "gene_pathway_membership_ensembl.csv",
  "gene_symbol_mapping_audit.csv",
  "gene_symbol_mapping_strict.csv",
  "core_gene_frequency_by_method.csv",
  "core_gene_frequency_overall.csv",
  "core_genes_with_differential_expression.csv",
  "top_20_core_genes.csv",
  "mapping_status_summary.csv",
  "DE_integration_status_summary.csv",
  "core_gene_summary.csv",
  "TRANSCRIPTOME_X_CORE_GENE_REPORT.md"
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

print(validation, row.names = FALSE)

if (!all(validation$exists)) {
  stop("Validazione output non superata.")
}

if (nrow(gene_membership) == 0) {
  stop("La tabella gene membership è vuota.")
}

if (nrow(integrated_core_genes) == 0) {
  stop("Nessun core gene integrato.")
}

# Specific check for the three previously unresolved genes.
target_ids <- c(
  "ENSG00000188522",
  "ENSG00000112697",
  "ENSG00000146540"
)

target_audit <- integrated_core_genes[
  integrated_core_genes$ensembl_id %in% target_ids,
  c(
    "comparison",
    "enrichment_symbol",
    "de_gene_name",
    "ensembl_id",
    "logFC",
    "FDR",
    "de_significant",
    "de_mapping_status"
  ),
  drop = FALSE
]

write.csv(
  target_audit,
  file.path(output_dir, "three_gene_mapping_recovery_audit.csv"),
  row.names = FALSE
)

cat("\n========== THREE GENE RECOVERY AUDIT ==========\n")
print(target_audit, row.names = FALSE)

cat("\n============================================================\n")
cat("CORE GENE ANALYSIS V2 COMPLETED\n")
cat("============================================================\n")

cat(
  "\nGeni distinti nei pathway:",
  length(unique(gene_membership$gene_key)),
  "\n"
)

cat(
  "Associazioni gene-term:",
  nrow(gene_membership),
  "\n"
)

cat(
  "Core genes integrati:",
  nrow(integrated_core_genes),
  "\n"
)

cat("\nCartella output:\n")
cat(output_dir, "\n")

cat("\nCORE GENE ANALYSIS V2 VALIDATED.\n")
