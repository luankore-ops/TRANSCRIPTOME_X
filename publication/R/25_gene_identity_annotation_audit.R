
# ============================================================
# TRANSCRIPTOME-X - GENE IDENTITY AND ANNOTATION AUDIT
# Script: 25_gene_identity_annotation_audit.R
#
# Obiettivi:
# 1. Verificare l'identita dei 533 geni concordanti di R24
# 2. Verificare associazioni Ensembl ID - gene symbol
# 3. Identificare simboli associati a piu Ensembl ID
# 4. Approfondire il caso CRLF2
# 5. Verificare la corrispondenza con org.Hs.eg.db
# 6. Verificare NRP1 nei due confronti DE
# 7. Validare FDR e direzione dei 533 geni
#
# Chiave primaria: Ensembl stable ID senza suffisso versione
#
# Input:
# - R24 all_genes_cross_knockdown.csv
# - R24 genes_concordant_both_knockdowns.csv
# - DE annotated originali
#
# Output indipendenti in:
# data/processed/differential_expression/
# cross_knockdown_annotation_audit/
#
# Non modifica script o output precedenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - GENE IDENTITY AND ANNOTATION AUDIT\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

project_root <- getwd()

de_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "annotated"
)

r24_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "cross_knockdown_concordance"
)

output_dir <- file.path(
  "data", "processed",
  "differential_expression",
  "cross_knockdown_annotation_audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

file_all <- file.path(
  r24_dir,
  "all_genes_cross_knockdown.csv"
)

file_concordant <- file.path(
  r24_dir,
  "genes_concordant_both_knockdowns.csv"
)

file_de_3 <- file.path(
  de_dir,
  "siNRP1_3_vs_siNT_annotated.csv"
)

file_de_5 <- file.path(
  de_dir,
  "siNRP1_5_vs_siNT_annotated.csv"
)

required_files <- c(
  file_all,
  file_concordant,
  file_de_3,
  file_de_5
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

# ------------------------------------------------------------
# 2. DIPENDENZE
# ------------------------------------------------------------

cat("\n1. VERIFICA DIPENDENZE\n")
cat("------------------------------------------------------------\n")

required_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Pacchetto richiesto non installato: ",
      pkg
    )
  }
  
  cat(
    pkg,
    ":",
    as.character(packageVersion(pkg)),
    "\n"
  )
}

library(AnnotationDbi)
library(org.Hs.eg.db)

annotation_db <- org.Hs.eg.db

if (!"ENSEMBL" %in% AnnotationDbi::keytypes(annotation_db)) {
  stop("ENSEMBL non disponibile come keytype in org.Hs.eg.db.")
}

cat(
  "Annotation database:",
  "org.Hs.eg.db",
  "\n"
)

# ------------------------------------------------------------
# 3. IMPORTAZIONE
# ------------------------------------------------------------

cat("\n2. IMPORTAZIONE DATASET\n")
cat("------------------------------------------------------------\n")

all_genes <- read.csv(
  file_all,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

concordant <- read.csv(
  file_concordant,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

de_3 <- read.csv(
  file_de_3,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

de_5 <- read.csv(
  file_de_5,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("R24 all genes:", nrow(all_genes), "\n")
cat("R24 concordant:", nrow(concordant), "\n")
cat("DE siNRP1_3:", nrow(de_3), "\n")
cat("DE siNRP1_5:", nrow(de_5), "\n")

# ------------------------------------------------------------
# 4. ISPEZIONE E VALIDAZIONE COLONNE
# ------------------------------------------------------------

cat("\n3. ISPEZIONE COLONNE REALI\n")
cat("------------------------------------------------------------\n")

cat("\nR24 all genes:\n")
print(names(all_genes))

cat("\nR24 concordant:\n")
print(names(concordant))

cat("\nDE siNRP1_3:\n")
print(names(de_3))

cat("\nDE siNRP1_5:\n")
print(names(de_5))

required_r24 <- c(
  "ensembl_id",
  "gene_name",
  "logFC_3",
  "FDR_3",
  "logFC_5",
  "FDR_5"
)

required_de <- c(
  "gene_id",
  "gene_name",
  "logFC",
  "FDR"
)

missing_r24 <- setdiff(
  required_r24,
  names(concordant)
)

if (length(missing_r24) > 0) {
  stop(
    "Colonne mancanti nel dataset concordant: ",
    paste(missing_r24, collapse = ", ")
  )
}

for (dataset_name in c("de_3", "de_5")) {
  
  dataset <- get(dataset_name)
  
  missing <- setdiff(
    required_de,
    names(dataset)
  )
  
  if (length(missing) > 0) {
    stop(
      "Colonne mancanti in ",
      dataset_name,
      ": ",
      paste(missing, collapse = ", ")
    )
  }
}

# ------------------------------------------------------------
# 5. NORMALIZZAZIONE ENSEMBL
# ------------------------------------------------------------

cat("\n4. NORMALIZZAZIONE ENSEMBL ID\n")
cat("------------------------------------------------------------\n")

normalize_ensembl <- function(x) {
  
  x <- trimws(as.character(x))
  
  sub(
    "\\.[0-9]+$",
    "",
    x
  )
}

concordant$ensembl_id <- normalize_ensembl(
  concordant$ensembl_id
)

all_genes$ensembl_id <- normalize_ensembl(
  all_genes$ensembl_id
)

de_3$ensembl_id <- normalize_ensembl(
  de_3$gene_id
)

de_5$ensembl_id <- normalize_ensembl(
  de_5$gene_id
)

if (anyNA(concordant$ensembl_id) ||
    any(!nzchar(concordant$ensembl_id))) {
  stop("Ensembl ID mancanti nel dataset concordant.")
}

if (anyDuplicated(concordant$ensembl_id) > 0) {
  stop("Ensembl ID duplicati nel dataset concordant.")
}

if (anyDuplicated(all_genes$ensembl_id) > 0) {
  stop("Ensembl ID duplicati nel dataset R24 all genes.")
}

if (anyDuplicated(de_3$ensembl_id) > 0 ||
    anyDuplicated(de_5$ensembl_id) > 0) {
  stop("Ensembl ID duplicati nei dataset DE originali.")
}

cat(
  "Ensembl univoci nei concordant:",
  length(unique(concordant$ensembl_id)),
  "\n"
)

# ------------------------------------------------------------
# 6. VALIDAZIONE APPARTENENZA AI 533
# ------------------------------------------------------------

cat("\n5. VALIDAZIONE APPARTENENZA AI 533 GENI\n")
cat("------------------------------------------------------------\n")

if (nrow(concordant) != 533) {
  warning(
    "Il numero di geni concordanti non corrisponde a 533: ",
    nrow(concordant)
  )
}

missing_from_all <- setdiff(
  concordant$ensembl_id,
  all_genes$ensembl_id
)

if (length(missing_from_all) > 0) {
  stop(
    "Alcuni geni concordanti non sono presenti in all_genes_cross_knockdown.csv."
  )
}

# ------------------------------------------------------------
# 7. VALIDAZIONE SIGNIFICATIVITA E DIREZIONE
# ------------------------------------------------------------

cat("\n6. VALIDAZIONE FDR E DIREZIONE\n")
cat("------------------------------------------------------------\n")

FDR_THRESHOLD <- 0.05

concordant$significant_3 <- (
  !is.na(concordant$FDR_3) &
    concordant$FDR_3 < FDR_THRESHOLD
)

concordant$significant_5 <- (
  !is.na(concordant$FDR_5) &
    concordant$FDR_5 < FDR_THRESHOLD
)

concordant$direction_3 <- ifelse(
  concordant$logFC_3 > 0,
  "UP",
  ifelse(concordant$logFC_3 < 0, "DOWN", "ZERO")
)

concordant$direction_5 <- ifelse(
  concordant$logFC_5 > 0,
  "UP",
  ifelse(concordant$logFC_5 < 0, "DOWN", "ZERO")
)

concordant$direction_consistent <- (
  concordant$direction_3 == concordant$direction_5 &
    concordant$direction_3 != "ZERO"
)

concordant$DE_validation <- (
  concordant$significant_3 &
    concordant$significant_5 &
    concordant$direction_consistent
)

cat(
  "Significativi siNRP1_3:",
  sum(concordant$significant_3),
  "\n"
)

cat(
  "Significativi siNRP1_5:",
  sum(concordant$significant_5),
  "\n"
)

cat(
  "Direzione coerente:",
  sum(concordant$direction_consistent),
  "\n"
)

cat(
  "Validati complessivamente:",
  sum(concordant$DE_validation),
  "\n"
)

if (!all(concordant$DE_validation)) {
  warning(
    "Alcuni geni non superano la validazione indipendente dei campi DE."
  )
}

write.csv(
  concordant,
  file.path(
    output_dir,
    "R25_concordant_genes_DE_validation.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. ANNOTAZIONE CON org.Hs.eg.db
# ------------------------------------------------------------

cat("\n7. ANNOTAZIONE org.Hs.eg.db\n")
cat("------------------------------------------------------------\n")

ensembl_keys <- unique(
  concordant$ensembl_id
)

annotation_raw <- AnnotationDbi::select(
  annotation_db,
  keys = ensembl_keys,
  keytype = "ENSEMBL",
  columns = c(
    "SYMBOL",
    "ENTREZID",
    "GENENAME"
  )
)

annotation_raw <- annotation_raw[
  !is.na(annotation_raw$ENSEMBL),
  ,
  drop = FALSE
]

annotation_raw <- unique(annotation_raw)

cat(
  "Righe restituite da AnnotationDbi:",
  nrow(annotation_raw),
  "\n"
)

# Funzione per aggregare valori distinti senza perdere
# le associazioni multiple.

collapse_unique <- function(x) {
  
  x <- unique(
    trimws(as.character(x))
  )
  
  x <- x[
    !is.na(x) &
      nzchar(x)
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    sort(x),
    collapse = ";"
  )
}

annotation_summary <- do.call(
  rbind,
  lapply(ensembl_keys, function(id) {
    
    current <- annotation_raw[
      annotation_raw$ENSEMBL == id,
      ,
      drop = FALSE
    ]
    
    data.frame(
      ensembl_id = id,
      annotation_symbols = collapse_unique(current$SYMBOL),
      annotation_entrez = collapse_unique(current$ENTREZID),
      annotation_genenames = collapse_unique(current$GENENAME),
      annotation_rows = nrow(current),
      stringsAsFactors = FALSE
    )
  })
)

annotation_summary$symbol_count <- vapply(
  strsplit(
    ifelse(
      is.na(annotation_summary$annotation_symbols),
      "",
      annotation_summary$annotation_symbols
    ),
    ";",
    fixed = TRUE
  ),
  function(x) {
    sum(nzchar(x))
  },
  integer(1)
)

annotation_summary$annotation_status <- ifelse(
  annotation_summary$symbol_count == 0,
  "NOT_MAPPED",
  ifelse(
    annotation_summary$symbol_count > 1,
    "MULTIPLE_SYMBOLS",
    "SINGLE_SYMBOL"
  )
)

# ------------------------------------------------------------
# 9. CONFRONTO SIMBOLO RNA-SEQ VS ANNOTAZIONE
# ------------------------------------------------------------

cat("\n8. CONFRONTO GENE SYMBOL\n")
cat("------------------------------------------------------------\n")

concordant$gene_name <- trimws(
  as.character(concordant$gene_name)
)

annotation_audit <- merge(
  concordant,
  annotation_summary,
  by = "ensembl_id",
  all.x = TRUE,
  sort = FALSE
)

annotation_audit$symbol_match_status <- vapply(
  seq_len(nrow(annotation_audit)),
  function(i) {
    
    source_symbol <- annotation_audit$gene_name[i]
    mapped_symbols <- annotation_audit$annotation_symbols[i]
    
    if (is.na(source_symbol) || !nzchar(source_symbol)) {
      return("SOURCE_SYMBOL_MISSING")
    }
    
    if (is.na(mapped_symbols) || !nzchar(mapped_symbols)) {
      return("ANNOTATION_NOT_AVAILABLE")
    }
    
    symbols <- strsplit(
      mapped_symbols,
      ";",
      fixed = TRUE
    )[[1]]
    
    if (source_symbol %in% symbols) {
      return("MATCH")
    }
    
    if (length(symbols) > 1) {
      return("MULTIPLE_SYMBOLS_NO_MATCH")
    }
    
    "MISMATCH"
  },
  character(1)
)

annotation_status_summary <- as.data.frame(
  table(annotation_audit$annotation_status)
)

names(annotation_status_summary) <- c(
  "annotation_status",
  "gene_count"
)

symbol_match_summary <- as.data.frame(
  table(annotation_audit$symbol_match_status)
)

names(symbol_match_summary) <- c(
  "symbol_match_status",
  "gene_count"
)

cat("\nStato annotazione:\n")
print(annotation_status_summary, row.names = FALSE)

cat("\nConfronto simboli:\n")
print(symbol_match_summary, row.names = FALSE)

write.csv(
  annotation_audit,
  file.path(
    output_dir,
    "R25_gene_annotation_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  annotation_status_summary,
  file.path(
    output_dir,
    "R25_annotation_status_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  symbol_match_summary,
  file.path(
    output_dir,
    "R25_symbol_match_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. SIMBOLI ASSOCIATI A PIU ENSEMBL ID
# ------------------------------------------------------------

cat("\n9. SIMBOLI ASSOCIATI A PIU ENSEMBL ID\n")
cat("------------------------------------------------------------\n")

symbol_groups <- split(
  annotation_audit$ensembl_id,
  annotation_audit$gene_name
)

duplicate_source_symbols <- names(
  Filter(
    function(x) length(unique(x)) > 1,
    symbol_groups
  )
)

duplicate_symbol_report <- do.call(
  rbind,
  lapply(duplicate_source_symbols, function(symbol) {
    
    ids <- sort(
      unique(
        symbol_groups[[symbol]]
      )
    )
    
    current <- annotation_audit[
      annotation_audit$gene_name == symbol,
      ,
      drop = FALSE
    ]
    
    data.frame(
      gene_name = symbol,
      ensembl_count = length(ids),
      ensembl_ids = paste(ids, collapse = ";"),
      directions = paste(
        unique(current$direction_3),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  })
)

if (is.null(duplicate_symbol_report)) {
  
  duplicate_symbol_report <- data.frame(
    gene_name = character(),
    ensembl_count = integer(),
    ensembl_ids = character(),
    directions = character()
  )
}

duplicate_symbol_report <- duplicate_symbol_report[
  order(
    -duplicate_symbol_report$ensembl_count,
    duplicate_symbol_report$gene_name
  ),
  ,
  drop = FALSE
]

cat(
  "Simboli associati a piu Ensembl ID:",
  nrow(duplicate_symbol_report),
  "\n"
)

if (nrow(duplicate_symbol_report) > 0) {
  print(
    duplicate_symbol_report,
    row.names = FALSE
  )
}

write.csv(
  duplicate_symbol_report,
  file.path(
    output_dir,
    "R25_duplicate_gene_symbols.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. CONTROLLO SPECIFICO CRLF2
# ------------------------------------------------------------

cat("\n10. APPROFONDIMENTO CRLF2\n")
cat("------------------------------------------------------------\n")

crlf2_ids <- c(
  "ENSG00000205755",
  "ENSG00000292363"
)

crlf2_report <- annotation_audit[
  annotation_audit$ensembl_id %in% crlf2_ids |
    annotation_audit$gene_name == "CRLF2" |
    annotation_audit$annotation_symbols == "CRLF2",
  ,
  drop = FALSE
]

cat(
  "Record CRLF2 nei 533 concordanti:",
  nrow(crlf2_report),
  "\n"
)

if (nrow(crlf2_report) > 0) {
  
  print(
    crlf2_report[
      ,
      intersect(
        c(
          "ensembl_id",
          "gene_name",
          "annotation_symbols",
          "annotation_entrez",
          "logFC_3",
          "FDR_3",
          "logFC_5",
          "FDR_5",
          "direction_3",
          "direction_5",
          "symbol_match_status"
        ),
        names(crlf2_report)
      )
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "Nessuno dei due identificativi CRLF2 risulta nei 533 concordanti.\n"
  )
}

write.csv(
  crlf2_report,
  file.path(
    output_dir,
    "R25_CRLF2_dedicated_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. CONTROLLO NRP1
# ------------------------------------------------------------

cat("\n11. CONTROLLO NRP1\n")
cat("------------------------------------------------------------\n")

NRP1_ID <- "ENSG00000099250"

extract_nrp1 <- function(dataset, comparison_name) {
  
  result <- dataset[
    dataset$ensembl_id == NRP1_ID,
    ,
    drop = FALSE
  ]
  
  if (nrow(result) == 0) {
    
    return(
      data.frame(
        comparison = comparison_name,
        ensembl_id = NRP1_ID,
        gene_name = NA_character_,
        logFC = NA_real_,
        FDR = NA_real_,
        significant = FALSE,
        direction = "NOT_FOUND",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(result) > 1) {
    stop("NRP1 associato a piu righe nel dataset ", comparison_name)
  }
  
  data.frame(
    comparison = comparison_name,
    ensembl_id = NRP1_ID,
    gene_name = result$gene_name,
    logFC = result$logFC,
    FDR = result$FDR,
    significant = (
      !is.na(result$FDR) &&
        result$FDR < FDR_THRESHOLD
    ),
    direction = ifelse(
      result$logFC > 0,
      "UP",
      ifelse(result$logFC < 0, "DOWN", "ZERO")
    ),
    stringsAsFactors = FALSE
  )
}

nrp1_report <- rbind(
  extract_nrp1(
    de_3,
    "siNRP1_3_vs_siNT"
  ),
  extract_nrp1(
    de_5,
    "siNRP1_5_vs_siNT"
  )
)

nrp1_report$included_in_533 <- (
  nrp1_report$ensembl_id %in%
    concordant$ensembl_id
)

cat("\nRisultati NRP1:\n")
print(nrp1_report, row.names = FALSE)

write.csv(
  nrp1_report,
  file.path(
    output_dir,
    "R25_NRP1_dedicated_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. IDENTIFICATIVI NON MAPPATI
# ------------------------------------------------------------

cat("\n12. IDENTIFICATIVI NON MAPPATI\n")
cat("------------------------------------------------------------\n")

unmapped <- annotation_audit[
  annotation_audit$annotation_status == "NOT_MAPPED" |
    is.na(annotation_audit$annotation_status),
  ,
  drop = FALSE
]

cat(
  "Geni senza annotazione SYMBOL:",
  nrow(unmapped),
  "\n"
)

write.csv(
  unmapped,
  file.path(
    output_dir,
    "R25_unmapped_ensembl_ids.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. METADATI ANNOTAZIONE
# ------------------------------------------------------------

annotation_metadata <- data.frame(
  package = c(
    "AnnotationDbi",
    "org.Hs.eg.db"
  ),
  version = c(
    as.character(packageVersion("AnnotationDbi")),
    as.character(packageVersion("org.Hs.eg.db"))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  annotation_metadata,
  file.path(
    output_dir,
    "R25_annotation_database_versions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. REPORT FINALE
# ------------------------------------------------------------

cat("\n13. GENERAZIONE REPORT\n")
cat("------------------------------------------------------------\n")

n_concordant <- nrow(concordant)

n_validated <- sum(
  concordant$DE_validation,
  na.rm = TRUE
)

n_match <- sum(
  annotation_audit$symbol_match_status == "MATCH",
  na.rm = TRUE
)

n_mismatch <- sum(
  annotation_audit$symbol_match_status == "MISMATCH",
  na.rm = TRUE
)

n_unmapped <- nrow(unmapped)

n_duplicate_symbols <- nrow(
  duplicate_symbol_report
)

nrp1_lines <- paste(
  paste0(
    "- ",
    nrp1_report$comparison,
    ": logFC = ",
    signif(nrp1_report$logFC, 5),
    "; FDR = ",
    signif(nrp1_report$FDR, 5),
    "; direzione = ",
    nrp1_report$direction
  ),
  collapse = "\n"
)

crlf2_lines <- if (nrow(crlf2_report) == 0) {
  "Nessuno dei due identificativi compare nei 533 concordanti."
} else {
  paste(
    paste0(
      "- ",
      crlf2_report$ensembl_id,
      " | ",
      crlf2_report$gene_name,
      " | ",
      crlf2_report$annotation_symbols
    ),
    collapse = "\n"
  )
}

report <- c(
  "# TRANSCRIPTOME-X - Gene Identity and Annotation Audit",
  "",
  "## Metodo",
  "",
  "- Chiave primaria: Ensembl stable ID",
  "- Rimozione suffisso versione Ensembl",
  "- Annotazione tramite AnnotationDbi e org.Hs.eg.db",
  "- Soglia DE: FDR < 0.05",
  "",
  "## Risultati principali",
  "",
  paste0("- Geni concordanti analizzati: ", n_concordant),
  paste0("- Geni che superano la validazione DE: ", n_validated),
  paste0("- Simboli corrispondenti all'annotazione: ", n_match),
  paste0("- Simboli con mismatch: ", n_mismatch),
  paste0("- Geni senza annotazione: ", n_unmapped),
  paste0("- Simboli associati a piu Ensembl ID: ", n_duplicate_symbols),
  "",
  "## Caso CRLF2",
  "",
  crlf2_lines,
  "",
  "## Controllo NRP1",
  "",
  nrp1_lines,
  "",
  "## Interpretazione",
  "",
  "Un simbolo genico associato a piu Ensembl ID non implica",
  "necessariamente un errore: possono esistere identificativi",
  "distinti che condividono un simbolo o rappresentano annotazioni",
  "diverse. Gli identificativi non vengono collassati.",
  "",
  "Le differenze tra gene_name e org.Hs.eg.db devono essere",
  "interpretate considerando le versioni e le fonti di annotazione.",
  "",
  "La corrispondenza con org.Hs.eg.db e un controllo di annotazione,",
  "non una validazione sperimentale dell'identita biologica.",
  "",
  "La validazione dei 533 geni utilizza gli stessi risultati DE",
  "di R24 e non costituisce una replica sperimentale."
)

writeLines(
  report,
  file.path(
    output_dir,
    "TRANSCRIPTOME_X_GENE_IDENTITY_ANNOTATION_AUDIT.md"
  )
)

# ------------------------------------------------------------
# 16. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

cat("\n14. VALIDAZIONE OUTPUT\n")
cat("------------------------------------------------------------\n")

output_files <- c(
  "R25_concordant_genes_DE_validation.csv",
  "R25_gene_annotation_audit.csv",
  "R25_annotation_status_summary.csv",
  "R25_symbol_match_summary.csv",
  "R25_duplicate_gene_symbols.csv",
  "R25_CRLF2_dedicated_audit.csv",
  "R25_NRP1_dedicated_audit.csv",
  "R25_unmapped_ensembl_ids.csv",
  "R25_annotation_database_versions.csv",
  "TRANSCRIPTOME_X_GENE_IDENTITY_ANNOTATION_AUDIT.md"
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
# 17. RIEPILOGO FINALE
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE IDENTITY AND ANNOTATION AUDIT COMPLETED\n")
cat("============================================================\n")

cat("\nRIEPILOGO\n")
cat("------------------------------------------------------------\n")

cat("Geni concordanti analizzati:", n_concordant, "\n")
cat("Validazione DE superata:", n_validated, "\n")
cat("Symbol match:", n_match, "\n")
cat("Symbol mismatch:", n_mismatch, "\n")
cat("Ensembl non mappati:", n_unmapped, "\n")
cat("Simboli con piu Ensembl ID:", n_duplicate_symbols, "\n")

cat("\nNRP1\n")
print(nrp1_report, row.names = FALSE)

cat("\nCRLF2\n")
if (nrow(crlf2_report) > 0) {
  print(
    crlf2_report[
      ,
      c(
        "ensembl_id",
        "gene_name",
        "annotation_symbols",
        "logFC_3",
        "FDR_3",
        "logFC_5",
        "FDR_5"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
} else {
  cat("Nessun record CRLF2 nei 533 concordanti.\n")
}

cat(
  "\nOutput directory:\n",
  normalizePath(output_dir, winslash = "/"),
  "\n"
)

cat("\nR25 AUDIT VALIDATED.\n")
