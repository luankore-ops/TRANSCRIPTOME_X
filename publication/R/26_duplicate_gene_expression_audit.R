
# ============================================================
# TRANSCRIPTOME-X - DUPLICATE GENE AND EXPRESSION AUDIT
# Script: 26_duplicate_gene_expression_audit.R
#
# Obiettivi:
# 1. Analizzare i conteggi RSEM originali dei 9 campioni
# 2. Verificare CRLF2, CSF2RA e PLCXD1
# 3. Identificare profili di espressione identici
# 4. Identificare simboli associati a piu Ensembl ID
# 5. Verificare la corrispondenza con i risultati DE
# 6. Produrre un report conclusivo
#
# Importante:
# - Non arrotonda i conteggi frazionari
# - Non elimina geni
# - Non modifica dati o output precedenti
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - DUPLICATE GENE EXPRESSION AUDIT\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

counts_file <- file.path(
  "data", "processed", "counts",
  "GSE266566_MDA_MB_231_selected_counts.tsv"
)

de_dir <- file.path(
  "data", "processed",
  "differential_expression", "annotated"
)

output_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "duplicate_gene_expression_audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

de_file_3 <- file.path(
  de_dir,
  "siNRP1_3_vs_siNT_annotated.csv"
)

de_file_5 <- file.path(
  de_dir,
  "siNRP1_5_vs_siNT_annotated.csv"
)

required_files <- c(
  counts_file,
  de_file_3,
  de_file_5
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "File mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

sample_ids <- c(
  "counts.MDAMB231.siNT.RNA.1",
  "counts.MDAMB231.siNT.RNA.2",
  "counts.MDAMB231.siNT.RNA.3",
  "counts.MDAMB231.siNRP1.3.RNA.1",
  "counts.MDAMB231.siNRP1.3.RNA.2",
  "counts.MDAMB231.siNRP1.3.RNA.3",
  "counts.MDAMB231.siNRP1.5.RNA.1",
  "counts.MDAMB231.siNRP1.5.RNA.2",
  "counts.MDAMB231.siNRP1.5.RNA.3"
)

target_symbols <- c(
  "CRLF2",
  "CSF2RA",
  "PLCXD1"
)

# ------------------------------------------------------------
# 2. IMPORTAZIONE CONTEGGI ORIGINALI
# ------------------------------------------------------------

cat("\n1. IMPORTAZIONE CONTEGGI RSEM\n")
cat("------------------------------------------------------------\n")

counts_data <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Righe originali:", nrow(counts_data), "\n")
cat("Colonne originali:", ncol(counts_data), "\n")

required_annotation <- c(
  "Ensembl.114.Gene.ID",
  "Ensembl.114.Gene.Name",
  "Ensembl.114.Gene.Type"
)

missing_annotation <- setdiff(
  required_annotation,
  names(counts_data)
)

if (length(missing_annotation) > 0) {
  stop(
    "Colonne di annotazione mancanti: ",
    paste(missing_annotation, collapse = ", ")
  )
}

missing_samples <- setdiff(
  sample_ids,
  names(counts_data)
)

if (length(missing_samples) > 0) {
  stop(
    "Campioni mancanti: ",
    paste(missing_samples, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. CREAZIONE MATRICE
# ------------------------------------------------------------

cat("\n2. CREAZIONE MATRICE DI ESPRESSIONE\n")
cat("------------------------------------------------------------\n")

gene_info <- data.frame(
  ensembl_original = as.character(
    counts_data$Ensembl.114.Gene.ID
  ),
  ensembl_id = sub(
    "\\.[0-9]+$",
    "",
    trimws(as.character(counts_data$Ensembl.114.Gene.ID))
  ),
  gene_name = trimws(
    as.character(counts_data$Ensembl.114.Gene.Name)
  ),
  gene_type = as.character(
    counts_data$Ensembl.114.Gene.Type
  ),
  stringsAsFactors = FALSE
)

counts <- as.matrix(
  counts_data[, sample_ids, drop = FALSE]
)

storage.mode(counts) <- "numeric"

if (anyNA(counts) || any(!is.finite(counts))) {
  stop("La matrice contiene conteggi mancanti o non finiti.")
}

if (any(counts < 0)) {
  stop("Sono presenti conteggi negativi.")
}

if (nrow(counts) != nrow(gene_info)) {
  stop("Dimensioni matrice e annotazioni non coerenti.")
}

cat("Geni:", nrow(counts), "\n")
cat("Campioni:", ncol(counts), "\n")
cat(
  "Conteggi frazionari:",
  sum(abs(counts - round(counts)) > 1e-8),
  "\n"
)

# ------------------------------------------------------------
# 4. VALIDAZIONE IDENTIFICATIVI
# ------------------------------------------------------------

cat("\n3. VALIDAZIONE IDENTIFICATIVI\n")
cat("------------------------------------------------------------\n")

cat(
  "Ensembl ID mancanti:",
  sum(
    is.na(gene_info$ensembl_id) |
      !nzchar(gene_info$ensembl_id)
  ),
  "\n"
)

duplicate_ensembl <- gene_info$ensembl_id[
  duplicated(gene_info$ensembl_id) |
    duplicated(gene_info$ensembl_id, fromLast = TRUE)
]

cat(
  "Ensembl ID duplicati dopo rimozione versione:",
  length(unique(duplicate_ensembl)),
  "\n"
)

# ------------------------------------------------------------
# 5. ANALISI PROFILI IDENTICI
# ------------------------------------------------------------

cat("\n4. IDENTIFICAZIONE PROFILI IDENTICI\n")
cat("------------------------------------------------------------\n")

# Firma costruita sui nove conteggi originali.
# sprintf conserva la precisione numerica disponibile.

profile_signature <- apply(
  counts,
  1,
  function(x) {
    paste(
      sprintf("%.17g", x),
      collapse = "|"
    )
  }
)

profile_groups <- split(
  seq_len(nrow(counts)),
  profile_signature
)

duplicate_groups <- Filter(
  function(x) length(x) > 1,
  profile_groups
)

cat(
  "Gruppi con profili identici:",
  length(duplicate_groups),
  "\n"
)

cat(
  "Geni coinvolti:",
  sum(lengths(duplicate_groups)),
  "\n"
)

if (length(duplicate_groups) > 0) {
  
  duplicate_profile_report <- do.call(
    rbind,
    lapply(seq_along(duplicate_groups), function(i) {
      
      idx <- duplicate_groups[[i]]
      
      data.frame(
        profile_group = paste0("PROFILE_", i),
        ensembl_ids = paste(
          gene_info$ensembl_id[idx],
          collapse = ";"
        ),
        gene_names = paste(
          unique(gene_info$gene_name[idx]),
          collapse = ";"
        ),
        gene_count = length(idx),
        total_counts = sum(counts[idx[1], ]),
        stringsAsFactors = FALSE
      )
    })
  )
  
} else {
  
  duplicate_profile_report <- data.frame(
    profile_group = character(),
    ensembl_ids = character(),
    gene_names = character(),
    gene_count = integer(),
    total_counts = numeric()
  )
}

write.csv(
  duplicate_profile_report,
  file.path(
    output_dir,
    "R26_identical_expression_profiles.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. AUDIT CRLF2, CSF2RA, PLCXD1
# ------------------------------------------------------------

cat("\n5. AUDIT GENI CON SIMBOLI MULTIPLI\n")
cat("------------------------------------------------------------\n")

candidate_indices <- which(
  gene_info$gene_name %in% target_symbols
)

candidate_genes <- gene_info[
  candidate_indices,
  ,
  drop = FALSE
]

candidate_counts <- counts[
  candidate_indices,
  ,
  drop = FALSE
]

candidate_report <- cbind(
  candidate_genes,
  as.data.frame(candidate_counts, check.names = FALSE)
)

candidate_report$profile_signature <- profile_signature[
  candidate_indices
]

candidate_report$profile_identical_within_symbol <- vapply(
  seq_along(candidate_indices),
  function(i) {
    
    same_symbol <- which(
      gene_info$gene_name ==
        candidate_genes$gene_name[i]
    )
    
    length(unique(
      profile_signature[same_symbol]
    )) == 1
  },
  logical(1)
)

write.csv(
  candidate_report,
  file.path(
    output_dir,
    "R26_candidate_gene_expression_profiles.csv"
  ),
  row.names = FALSE
)

for (symbol in target_symbols) {
  
  cat("\n", symbol, "\n", sep = "")
  
  idx <- which(
    gene_info$gene_name == symbol
  )
  
  if (length(idx) == 0) {
    cat("Non presente nella matrice originale.\n")
    next
  }
  
  print(
    data.frame(
      ensembl_id = gene_info$ensembl_id[idx],
      gene_name = gene_info$gene_name[idx],
      total_counts = rowSums(counts[idx, , drop = FALSE]),
      identical_profile = duplicated(
        profile_signature[idx]
      ) |
        duplicated(
          profile_signature[idx],
          fromLast = TRUE
        ),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
  if (length(idx) > 1) {
    
    cat(
      "Profili esattamente identici:",
      length(unique(profile_signature[idx])) == 1,
      "\n"
    )
    
    cat(
      "Differenza massima assoluta:",
      max(
        abs(
          counts[idx[1], ] -
            counts[idx[2], ]
        )
      ),
      "\n"
    )
  }
}

# ------------------------------------------------------------
# 7. SIMBOLI MULTIPLI NELLA MATRICE COMPLETA
# ------------------------------------------------------------

cat("\n6. SIMBOLI ASSOCIATI A PIU ENSEMBL ID\n")
cat("------------------------------------------------------------\n")

valid_symbols <- gene_info$gene_name[
  !is.na(gene_info$gene_name) &
    nzchar(gene_info$gene_name)
]

symbol_id_counts <- aggregate(
  gene_info$ensembl_id[
    !is.na(gene_info$gene_name) &
      nzchar(gene_info$gene_name)
  ],
  by = list(
    gene_name = valid_symbols
  ),
  FUN = function(x) length(unique(x))
)

names(symbol_id_counts)[2] <- "ensembl_count"

multiple_symbols <- symbol_id_counts[
  symbol_id_counts$ensembl_count > 1,
  ,
  drop = FALSE
]

multiple_symbols <- multiple_symbols[
  order(-multiple_symbols$ensembl_count),
  ,
  drop = FALSE
]

cat(
  "Simboli associati a piu Ensembl ID:",
  nrow(multiple_symbols),
  "\n"
)

write.csv(
  multiple_symbols,
  file.path(
    output_dir,
    "R26_multiple_symbols_full_matrix.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. CONFRONTO CON I RISULTATI DIFFERENTIAL EXPRESSION
# ------------------------------------------------------------

cat("\n7. CONFRONTO CON DIFFERENTIAL EXPRESSION\n")
cat("------------------------------------------------------------\n")

de_3 <- read.csv(
  de_file_3,
  stringsAsFactors = FALSE
)

de_5 <- read.csv(
  de_file_5,
  stringsAsFactors = FALSE
)

de_3$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  trimws(as.character(de_3$gene_id))
)

de_5$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  trimws(as.character(de_5$gene_id))
)

de_3_small <- de_3[
  ,
  c("ensembl_id", "gene_name", "logFC", "FDR")
]

names(de_3_small) <- c(
  "ensembl_id",
  "gene_name_3",
  "logFC_3",
  "FDR_3"
)

de_5_small <- de_5[
  ,
  c("ensembl_id", "gene_name", "logFC", "FDR")
]

names(de_5_small) <- c(
  "ensembl_id",
  "gene_name_5",
  "logFC_5",
  "FDR_5"
)

de_combined <- merge(
  de_3_small,
  de_5_small,
  by = "ensembl_id",
  all = TRUE
)

candidate_de <- merge(
  candidate_report,
  de_combined,
  by = "ensembl_id",
  all.x = TRUE
)

write.csv(
  candidate_de,
  file.path(
    output_dir,
    "R26_candidate_genes_DE_comparison.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. CONFRONTO DIRETTO DEI PROFILI
# ------------------------------------------------------------

cat("\n8. CONFRONTO PROFILI CANDIDATI\n")
cat("------------------------------------------------------------\n")

candidate_pair_results <- list()

for (symbol in target_symbols) {
  
  idx <- which(
    gene_info$gene_name == symbol
  )
  
  if (length(idx) < 2) {
    next
  }
  
  pairs <- combn(idx, 2)
  
  for (j in seq_len(ncol(pairs))) {
    
    i1 <- pairs[1, j]
    i2 <- pairs[2, j]
    
    x <- counts[i1, ]
    y <- counts[i2, ]
    
    candidate_pair_results[[length(candidate_pair_results) + 1]] <-
      data.frame(
        gene_name = symbol,
        ensembl_id_1 = gene_info$ensembl_id[i1],
        ensembl_id_2 = gene_info$ensembl_id[i2],
        identical_profile = identical(
          as.numeric(x),
          as.numeric(y)
        ),
        max_absolute_difference = max(abs(x - y)),
        pearson_correlation = if (
          sd(x) > 0 && sd(y) > 0
        ) {
          cor(x, y, method = "pearson")
        } else {
          NA_real_
        },
        total_count_1 = sum(x),
        total_count_2 = sum(y),
        stringsAsFactors = FALSE
      )
  }
}

if (length(candidate_pair_results) > 0) {
  
  candidate_pair_report <- do.call(
    rbind,
    candidate_pair_results
  )
  
} else {
  
  candidate_pair_report <- data.frame()
}

print(candidate_pair_report, row.names = FALSE)

write.csv(
  candidate_pair_report,
  file.path(
    output_dir,
    "R26_candidate_pairwise_comparison.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. CONTROLLI FINALI
# ------------------------------------------------------------

cat("\n9. CONTROLLI FINALI\n")
cat("------------------------------------------------------------\n")

candidate_summary <- do.call(
  rbind,
  lapply(target_symbols, function(symbol) {
    
    idx <- which(
      gene_info$gene_name == symbol
    )
    
    data.frame(
      gene_name = symbol,
      ensembl_count = length(idx),
      unique_profiles = length(
        unique(profile_signature[idx])
      ),
      identical_profiles = (
        length(idx) > 1 &&
          length(unique(profile_signature[idx])) == 1
      ),
      stringsAsFactors = FALSE
    )
  })
)

print(candidate_summary, row.names = FALSE)

write.csv(
  candidate_summary,
  file.path(
    output_dir,
    "R26_candidate_gene_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. REPORT
# ------------------------------------------------------------

cat("\n10. REPORT FINALE\n")
cat("------------------------------------------------------------\n")

n_duplicate_groups <- length(duplicate_groups)

n_duplicate_genes <- sum(
  lengths(duplicate_groups)
)

candidate_lines <- paste(
  paste0(
    "- ",
    candidate_summary$gene_name,
    ": ",
    candidate_summary$ensembl_count,
    " Ensembl ID; ",
    candidate_summary$unique_profiles,
    " profili distinti"
  ),
  collapse = "\n"
)

report <- c(
  "# TRANSCRIPTOME-X - Duplicate Gene Expression Audit",
  "",
  "## Dataset",
  "",
  "- Fonte: RSEM original selected counts",
  "- Campioni: 9",
  paste0("- Geni: ", nrow(counts)),
  "- Conteggi frazionari mantenuti senza arrotondamento",
  "",
  "## Profili identici",
  "",
  paste0("- Gruppi di profili identici: ", n_duplicate_groups),
  paste0("- Geni coinvolti: ", n_duplicate_genes),
  "",
  "## Geni con simboli multipli",
  "",
  candidate_lines,
  "",
  "## Interpretazione",
  "",
  "Profili identici sui nove campioni indicano che due righe",
  "hanno gli stessi valori di conteggio nella matrice analizzata.",
  "",
  "Questo risultato non dimostra da solo se la duplicazione",
  "derivi dall'annotazione, dalla quantificazione o dalla",
  "rappresentazione di feature genomiche distinte.",
  "",
  "Non sono stati rimossi o aggregati geni.",
  "",
  "La matrice analizzata contiene i nove campioni selezionati",
  "e non tutti i campioni del dataset GEO originale."
)

writeLines(
  report,
  file.path(
    output_dir,
    "TRANSCRIPTOME_X_DUPLICATE_GENE_EXPRESSION_AUDIT.md"
  )
)

# ------------------------------------------------------------
# 12. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n11. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "R26_identical_expression_profiles.csv",
  "R26_candidate_gene_expression_profiles.csv",
  "R26_multiple_symbols_full_matrix.csv",
  "R26_candidate_genes_DE_comparison.csv",
  "R26_candidate_pairwise_comparison.csv",
  "R26_candidate_gene_summary.csv",
  "TRANSCRIPTOME_X_DUPLICATE_GENE_EXPRESSION_AUDIT.md"
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
  stop("Validazione output fallita.")
}

# ------------------------------------------------------------
# 13. RIEPILOGO
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DUPLICATE GENE EXPRESSION AUDIT COMPLETED\n")
cat("============================================================\n")

cat("Geni analizzati:", nrow(counts), "\n")
cat("Campioni analizzati:", ncol(counts), "\n")
cat("Gruppi con profili identici:", n_duplicate_groups, "\n")
cat("Geni in profili identici:", n_duplicate_genes, "\n")

cat("\nCANDIDATI\n")
print(candidate_summary, row.names = FALSE)

cat(
  "\nOutput:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nR26 AUDIT COMPLETED.\n")
