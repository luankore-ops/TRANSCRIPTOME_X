
# ============================================================
# TRANSCRIPTOME-X
# Script: 32_functional_enrichment_interpretation.R
#
# FUNCTIONAL ENRICHMENT INTERPRETATION
# AND PATHWAY CONSOLIDATION
#
# Obiettivi:
# 1. Importare i risultati di arricchimento dello Script 31.
# 2. Identificare i termini significativi (BH-adjusted p < 0.05).
# 3. Consolidare termini ridondanti tramite overlap dei geni.
# 4. Confrontare i gruppi robusti e dipendenti dalla normalizzazione.
# 5. Confrontare i due contrasti sperimentali.
# 6. Produrre tabelle interpretabili e un report riepilogativo.
#
# Nota metodologica:
# La consolidazione è basata sulla sovrapposizione dei geni
# associati ai termini, NON sulla similarità semantica GO.
# I risultati originali non vengono modificati.
# ============================================================


# ------------------------------------------------------------
# 0. CONFIGURAZIONE
# ------------------------------------------------------------

PROJECT_DIR <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

INPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "discordant_gene_functional_analysis"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "functional_enrichment_interpretation"
)

FDR_THRESHOLD <- 0.05

# Due termini vengono collegati se:
# - Jaccard gene overlap >= 0.50
# oppure
# - almeno l'80% dei geni del termine più piccolo
# è condiviso con l'altro termine.
JACCARD_THRESHOLD <- 0.50
CONTAINMENT_THRESHOLD <- 0.80

CONTRASTS <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

GROUP_ORDER <- c(
  "COMMON_ALL_METHODS",
  "DISCORDANT_ANY_METHOD",
  "ALTERNATIVE_ONLY",
  "MIXED_METHOD_DISCORDANCE",
  "TMM_ONLY",
  "EXCLUSIVE_TMMwsp",
  "EXCLUSIVE_RLE",
  "EXCLUSIVE_UPPERQUARTILE"
)

ONTOLOGY_ORDER <- c(
  "BP",
  "CC",
  "MF",
  "Reactome"
)


# ------------------------------------------------------------
# 1. INIZIALIZZAZIONE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("FUNCTIONAL ENRICHMENT INTERPRETATION\n")
cat("============================================================\n\n")

if (!dir.exists(INPUT_DIR)) {
  stop(
    "Directory di input non trovata: ",
    INPUT_DIR
  )
}

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")
AUDIT_DIR <- file.path(OUTPUT_DIR, "audit")

for (directory in c(TABLE_DIR, REPORT_DIR, AUDIT_DIR)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("Directory input:\n", INPUT_DIR, "\n\n")
cat("Directory output:\n", OUTPUT_DIR, "\n\n")


# ------------------------------------------------------------
# 2. FUNZIONI DI SUPPORTO
# ------------------------------------------------------------

write_csv_safe <- function(x, filename) {
  
  if (is.null(x)) {
    return(invisible(NULL))
  }
  
  write.csv(
    x,
    file = filename,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}


first_existing_column <- function(data, candidates) {
  
  available <- candidates[candidates %in% names(data)]
  
  if (length(available) == 0) {
    return(NA_character_)
  }
  
  available[[1]]
}


clean_text <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[x == ""] <- NA_character_
  
  x
}


normalize_ontology <- function(x) {
  
  x <- tolower(trimws(as.character(x)))
  
  result <- rep(NA_character_, length(x))
  
  result[x %in% c("bp", "biological_process")] <- "BP"
  result[x %in% c("cc", "cellular_component")] <- "CC"
  result[x %in% c("mf", "molecular_function")] <- "MF"
  
  result[
    grepl("reactome", x, ignore.case = TRUE)
  ] <- "Reactome"
  
  result
}


normalize_group <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[x == ""] <- NA_character_
  
  x
}


parse_filename_metadata <- function(filename) {
  
  basename_file <- basename(filename)
  
  contrast <- NA_character_
  group <- NA_character_
  ontology <- NA_character_
  
  for (candidate in CONTRASTS) {
    
    if (grepl(candidate, basename_file, fixed = TRUE)) {
      contrast <- candidate
      break
    }
  }
  
  for (candidate in GROUP_ORDER) {
    
    if (grepl(candidate, basename_file, fixed = TRUE)) {
      group <- candidate
      break
    }
  }
  
  ontology_candidates <- c(
    "BP",
    "CC",
    "MF",
    "Reactome"
  )
  
  for (candidate in ontology_candidates) {
    
    if (grepl(candidate, basename_file, ignore.case = TRUE)) {
      ontology <- candidate
      break
    }
  }
  
  list(
    contrast = contrast,
    group = group,
    ontology = ontology
  )
}


parse_gene_ids <- function(x) {
  
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }
  
  ids <- unlist(
    strsplit(
      as.character(x),
      split = "[/;,|]+"
    ),
    use.names = FALSE
  )
  
  ids <- trimws(ids)
  ids <- ids[nzchar(ids)]
  ids <- unique(ids)
  
  ids
}


safe_numeric <- function(x) {
  
  suppressWarnings(
    as.numeric(as.character(x))
  )
}


safe_integer <- function(x) {
  
  result <- suppressWarnings(
    as.integer(as.character(x))
  )
  
  result
}


# ------------------------------------------------------------
# 3. IMPORTAZIONE DEI RISULTATI DI ENRICHMENT
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("IMPORTAZIONE RISULTATI\n")
cat("------------------------------------------------------------\n\n")

csv_files <- list.files(
  INPUT_DIR,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(csv_files) == 0) {
  stop("Nessun file CSV trovato nella directory di input.")
}

required_columns <- c(
  "ID",
  "Description",
  "geneID",
  "p.adjust"
)

imported_results <- list()
import_log <- list()

for (i in seq_along(csv_files)) {
  
  current_file <- csv_files[[i]]
  
  imported <- tryCatch(
    read.csv(
      current_file,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8"
    ),
    error = function(e) NULL
  )
  
  if (is.null(imported)) {
    
    import_log[[length(import_log) + 1]] <- data.frame(
      file = current_file,
      status = "READ_ERROR",
      rows = NA_integer_,
      reason = "Impossibile leggere il CSV",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  if (!all(required_columns %in% names(imported))) {
    
    import_log[[length(import_log) + 1]] <- data.frame(
      file = current_file,
      status = "SKIPPED_NON_ENRICHMENT",
      rows = nrow(imported),
      reason = "Colonne di arricchimento non riconosciute",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  if (nrow(imported) == 0) {
    
    import_log[[length(import_log) + 1]] <- data.frame(
      file = current_file,
      status = "EMPTY",
      rows = 0L,
      reason = "File senza risultati",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  metadata <- parse_filename_metadata(current_file)
  
  contrast_column <- first_existing_column(
    imported,
    c("contrast", "Contrast")
  )
  
  group_column <- first_existing_column(
    imported,
    c("group", "Group")
  )
  
  ontology_column <- first_existing_column(
    imported,
    c("ontology", "Ontology", "ONTOLOGY")
  )
  
  if (!is.na(contrast_column)) {
    imported$contrast <- clean_text(imported[[contrast_column]])
  } else {
    imported$contrast <- metadata$contrast
  }
  
  if (!is.na(group_column)) {
    imported$group <- normalize_group(imported[[group_column]])
  } else {
    imported$group <- metadata$group
  }
  
  if (!is.na(ontology_column)) {
    imported$ontology <- normalize_ontology(imported[[ontology_column]])
  } else {
    imported$ontology <- normalize_ontology(metadata$ontology)
  }
  
  imported$contrast[is.na(imported$contrast)] <- metadata$contrast
  imported$group[is.na(imported$group)] <- metadata$group
  imported$ontology[is.na(imported$ontology)] <- normalize_ontology(
    metadata$ontology
  )
  
  imported$source_file <- normalizePath(
    current_file,
    winslash = "/",
    mustWork = FALSE
  )
  
  imported$p.adjust <- safe_numeric(imported$p.adjust)
  
  if ("pvalue" %in% names(imported)) {
    imported$pvalue <- safe_numeric(imported$pvalue)
  } else {
    imported$pvalue <- NA_real_
  }
  
  if ("Count" %in% names(imported)) {
    imported$Count <- safe_integer(imported$Count)
  } else {
    imported$Count <- NA_integer_
  }
  
  imported$geneID <- as.character(imported$geneID)
  imported$ID <- as.character(imported$ID)
  imported$Description <- as.character(imported$Description)
  
  imported_results[[length(imported_results) + 1]] <- imported
  
  import_log[[length(import_log) + 1]] <- data.frame(
    file = current_file,
    status = "IMPORTED",
    rows = nrow(imported),
    reason = "",
    stringsAsFactors = FALSE
  )
}

if (length(import_log) > 0) {
  
  import_log_df <- do.call(
    rbind,
    import_log
  )
  
  write_csv_safe(
    import_log_df,
    file.path(AUDIT_DIR, "import_audit.csv")
  )
  
} else {
  
  import_log_df <- data.frame()
}

if (length(imported_results) == 0) {
  stop(
    "Nessun file di arricchimento valido trovato. ",
    "Verificare la struttura dei CSV dello Script 31."
  )
}

enrichment <- do.call(
  rbind,
  imported_results
)

rownames(enrichment) <- NULL

# Mantiene solo record con identificativi e metadati utilizzabili.
valid_rows <- (
  !is.na(enrichment$ID) &
    !is.na(enrichment$Description) &
    !is.na(enrichment$contrast) &
    !is.na(enrichment$group) &
    !is.na(enrichment$ontology)
)

enrichment <- enrichment[valid_rows, , drop = FALSE]

if (nrow(enrichment) == 0) {
  stop(
    "I CSV sono stati letti, ma nessuna riga contiene ",
    "tutti i metadati necessari."
  )
}

enrichment$gene_list <- lapply(
  enrichment$geneID,
  parse_gene_ids
)

cat("File CSV analizzati:", length(csv_files), "\n")
cat("File di arricchimento validi:", length(imported_results), "\n")
cat("Record importati:", nrow(enrichment), "\n")
cat(
  "Contrasti:",
  paste(sort(unique(enrichment$contrast)), collapse = ", "),
  "\n"
)
cat(
  "Gruppi:",
  length(unique(enrichment$group)),
  "\n"
)
cat(
  "Ontologie:",
  paste(sort(unique(enrichment$ontology)), collapse = ", "),
  "\n\n"
)


# ------------------------------------------------------------
# 4. INVENTARIO COMPLETO DEI TERMINI
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("INVENTARIO DEI TERMINI\n")
cat("------------------------------------------------------------\n\n")

enrichment$significant <- (
  !is.na(enrichment$p.adjust) &
    enrichment$p.adjust < FDR_THRESHOLD
)

enrichment$term_gene_count <- vapply(
  enrichment$gene_list,
  length,
  integer(1)
)

term_inventory <- enrichment[
  ,
  c(
    "contrast",
    "group",
    "ontology",
    "ID",
    "Description",
    "pvalue",
    "p.adjust",
    "qvalue",
    "Count",
    "term_gene_count",
    "significant",
    "source_file"
  )[c(
    "contrast",
    "group",
    "ontology",
    "ID",
    "Description",
    "pvalue",
    "p.adjust",
    "qvalue",
    "Count",
    "term_gene_count",
    "significant",
    "source_file"
  ) %in% names(enrichment)],
  drop = FALSE
]

term_inventory <- term_inventory[
  order(
    term_inventory$contrast,
    term_inventory$group,
    term_inventory$ontology,
    term_inventory$p.adjust,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write_csv_safe(
  term_inventory,
  file.path(TABLE_DIR, "complete_enrichment_inventory.csv")
)

significant_terms <- enrichment[
  enrichment$significant,
  ,
  drop = FALSE
]

write_csv_safe(
  significant_terms[
    ,
    setdiff(names(significant_terms), "gene_list"),
    drop = FALSE
  ],
  file.path(TABLE_DIR, "significant_enrichment_terms.csv")
)

cat("Termini complessivi:", nrow(enrichment), "\n")
cat("Termini significativi FDR < ", FDR_THRESHOLD, ": ",
    nrow(significant_terms), "\n\n", sep = "")


# ------------------------------------------------------------
# 5. RIEPILOGO PER CONTRASTO, GRUPPO E ONTOLOGIA
# ------------------------------------------------------------

summary_keys <- unique(
  enrichment[
    ,
    c("contrast", "group", "ontology"),
    drop = FALSE
  ]
)

summary_list <- vector(
  "list",
  nrow(summary_keys)
)

for (i in seq_len(nrow(summary_keys))) {
  
  current <- enrichment[
    enrichment$contrast == summary_keys$contrast[i] &
      enrichment$group == summary_keys$group[i] &
      enrichment$ontology == summary_keys$ontology[i],
    ,
    drop = FALSE
  ]
  
  summary_list[[i]] <- data.frame(
    contrast = summary_keys$contrast[i],
    group = summary_keys$group[i],
    ontology = summary_keys$ontology[i],
    total_terms = nrow(current),
    significant_terms = sum(current$significant, na.rm = TRUE),
    minimum_FDR = if (all(is.na(current$p.adjust))) {
      NA_real_
    } else {
      min(current$p.adjust, na.rm = TRUE)
    },
    stringsAsFactors = FALSE
  )
}

enrichment_summary <- do.call(
  rbind,
  summary_list
)

write_csv_safe(
  enrichment_summary,
  file.path(TABLE_DIR, "enrichment_summary_by_group.csv")
)


# ------------------------------------------------------------
# 6. CONSOLIDAMENTO DEI TERMINI RIDONDANTI
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CONSOLIDAMENTO DEI TERMINI SIGNIFICATIVI\n")
cat("------------------------------------------------------------\n\n")


# Union-find per identificare componenti connesse.
find_root <- function(parent, x) {
  
  while (parent[x] != x) {
    x <- parent[x]
  }
  
  x
}


union_sets <- function(parent, a, b) {
  
  root_a <- find_root(parent, a)
  root_b <- find_root(parent, b)
  
  if (root_a != root_b) {
    parent[root_b] <- root_a
  }
  
  parent
}


consolidate_terms <- function(data) {
  
  if (nrow(data) == 0) {
    return(data.frame())
  }
  
  data <- data[
    order(
      data$p.adjust,
      -ifelse(is.na(data$Count), 0, data$Count),
      data$ID,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  rownames(data) <- NULL
  
  n <- nrow(data)
  
  parent <- seq_len(n)
  
  gene_sets <- data$gene_list
  
  if (n >= 2) {
    
    for (i in seq_len(n - 1)) {
      
      genes_i <- gene_sets[[i]]
      
      if (length(genes_i) == 0) {
        next
      }
      
      for (j in seq.int(i + 1, n)) {
        
        genes_j <- gene_sets[[j]]
        
        if (length(genes_j) == 0) {
          next
        }
        
        intersection_size <- length(
          intersect(genes_i, genes_j)
        )
        
        if (intersection_size == 0) {
          next
        }
        
        union_size <- length(
          union(genes_i, genes_j)
        )
        
        jaccard <- intersection_size / union_size
        
        containment <- intersection_size / min(
          length(genes_i),
          length(genes_j)
        )
        
        if (
          jaccard >= JACCARD_THRESHOLD ||
          containment >= CONTAINMENT_THRESHOLD
        ) {
          
          parent <- union_sets(
            parent,
            i,
            j
          )
        }
      }
    }
  }
  
  roots <- vapply(
    seq_len(n),
    function(x) find_root(parent, x),
    integer(1)
  )
  
  cluster_ids <- match(
    roots,
    unique(roots)
  )
  
  data$cluster_id <- paste0(
    "CLUSTER_",
    sprintf("%03d", cluster_ids)
  )
  
  cluster_sizes <- table(cluster_ids)
  
  data$cluster_size <- as.integer(
    cluster_sizes[cluster_ids]
  )
  
  data$is_representative <- FALSE
  
  for (cluster in unique(cluster_ids)) {
    
    indices <- which(cluster_ids == cluster)
    
    # Il primo elemento è il rappresentante:
    # FDR più basso, poi Count più alto.
    representative <- indices[1]
    
    data$is_representative[representative] <- TRUE
  }
  
  data$cluster_member_IDs <- vapply(
    cluster_ids,
    function(cluster) {
      paste(
        data$ID[cluster_ids == cluster],
        collapse = " | "
      )
    },
    character(1)
  )
  
  data$cluster_member_descriptions <- vapply(
    cluster_ids,
    function(cluster) {
      paste(
        data$Description[cluster_ids == cluster],
        collapse = " | "
      )
    },
    character(1)
  )
  
  data$cluster_member_FDR <- vapply(
    cluster_ids,
    function(cluster) {
      paste(
        signif(data$p.adjust[cluster_ids == cluster], 4),
        collapse = " | "
      )
    },
    character(1)
  )
  
  data
}


consolidation_keys <- unique(
  significant_terms[
    ,
    c("contrast", "group", "ontology"),
    drop = FALSE
  ]
)

consolidated_list <- list()

if (nrow(significant_terms) > 0) {
  
  for (i in seq_len(nrow(consolidation_keys))) {
    
    current <- significant_terms[
      significant_terms$contrast == consolidation_keys$contrast[i] &
        significant_terms$group == consolidation_keys$group[i] &
        significant_terms$ontology == consolidation_keys$ontology[i],
      ,
      drop = FALSE
    ]
    
    consolidated <- consolidate_terms(current)
    
    if (nrow(consolidated) > 0) {
      consolidated_list[[length(consolidated_list) + 1]] <- consolidated
    }
  }
}

if (length(consolidated_list) > 0) {
  
  consolidated_terms <- do.call(
    rbind,
    consolidated_list
  )
  
  rownames(consolidated_terms) <- NULL
  
} else {
  
  consolidated_terms <- data.frame()
}

if (nrow(consolidated_terms) > 0) {
  
  consolidated_export <- consolidated_terms[
    ,
    setdiff(names(consolidated_terms), "gene_list"),
    drop = FALSE
  ]
  
  write_csv_safe(
    consolidated_export,
    file.path(TABLE_DIR, "redundancy_consolidated_terms.csv")
  )
  
  representative_terms <- consolidated_terms[
    consolidated_terms$is_representative,
    ,
    drop = FALSE
  ]
  
  write_csv_safe(
    representative_terms[
      ,
      setdiff(names(representative_terms), "gene_list"),
      drop = FALSE
    ],
    file.path(TABLE_DIR, "representative_functional_terms.csv")
  )
  
} else {
  
  representative_terms <- data.frame()
  
  write_csv_safe(
    data.frame(),
    file.path(TABLE_DIR, "redundancy_consolidated_terms.csv")
  )
  
  write_csv_safe(
    data.frame(),
    file.path(TABLE_DIR, "representative_functional_terms.csv")
  )
}


# ------------------------------------------------------------
# 7. CONFRONTO TRA GRUPPI DI ROBUSTEZZA
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CONFRONTO TRA GRUPPI DI ROBUSTEZZA\n")
cat("------------------------------------------------------------\n\n")


group_role <- function(group) {
  
  if (group == "COMMON_ALL_METHODS") {
    return("ROBUST_ALL_NORMALIZATIONS")
  }
  
  if (group == "DISCORDANT_ANY_METHOD") {
    return("NORMALIZATION_SENSITIVE")
  }
  
  if (group %in% c(
    "ALTERNATIVE_ONLY",
    "MIXED_METHOD_DISCORDANCE"
  )) {
    return("TMM_NOT_UNIFORMLY_SUPPORTED")
  }
  
  if (group == "TMM_ONLY") {
    return("TMM_SPECIFIC")
  }
  
  if (grepl("^EXCLUSIVE_", group)) {
    return("ALTERNATIVE_METHOD_SPECIFIC")
  }
  
  "OTHER"
}


if (nrow(enrichment) > 0) {
  
  enrichment$group_role <- vapply(
    enrichment$group,
    group_role,
    character(1)
  )
}

group_role_summary <- aggregate(
  significant ~ contrast + group_role + ontology,
  data = enrichment,
  FUN = sum,
  na.rm = TRUE
)

names(group_role_summary)[
  names(group_role_summary) == "significant"
] <- "significant_terms"

write_csv_safe(
  group_role_summary,
  file.path(TABLE_DIR, "significant_terms_by_robustness_role.csv")
)


# ------------------------------------------------------------
# 8. CONFRONTO DEI TERMINI TRA GRUPPI
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CONFRONTO DEI TERMINI TRA GRUPPI\n")
cat("------------------------------------------------------------\n\n")

if (nrow(significant_terms) > 0) {
  
  term_presence <- unique(
    significant_terms[
      ,
      c(
        "contrast",
        "ontology",
        "ID",
        "Description",
        "group"
      ),
      drop = FALSE
    ]
  )
  
  term_presence$group_role <- vapply(
    term_presence$group,
    group_role,
    character(1)
  )
  
  term_presence_summary <- aggregate(
    group ~ contrast + ontology + ID + Description,
    data = term_presence,
    FUN = function(x) paste(
      sort(unique(x)),
      collapse = " | "
    )
  )
  
  names(term_presence_summary)[
    names(term_presence_summary) == "group"
  ] <- "significant_groups"
  
  term_presence_summary$number_of_groups <- vapply(
    strsplit(
      term_presence_summary$significant_groups,
      " \\| "
    ),
    length,
    integer(1)
  )
  
  write_csv_safe(
    term_presence_summary,
    file.path(TABLE_DIR, "term_presence_across_groups.csv")
  )
  
  role_presence <- aggregate(
    group_role ~ contrast + ontology + ID + Description,
    data = term_presence,
    FUN = function(x) paste(
      sort(unique(x)),
      collapse = " | "
    )
  )
  
  names(role_presence)[
    names(role_presence) == "group_role"
  ] <- "significant_roles"
  
  write_csv_safe(
    role_presence,
    file.path(TABLE_DIR, "term_presence_across_robustness_roles.csv")
  )
}


# ------------------------------------------------------------
# 9. CONFRONTO TRA I DUE KNOCKDOWN
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CONFRONTO TRA CONTRASTI\n")
cat("------------------------------------------------------------\n\n")

if (nrow(significant_terms) > 0) {
  
  contrast_presence <- unique(
    significant_terms[
      ,
      c(
        "contrast",
        "ontology",
        "ID",
        "Description",
        "group"
      ),
      drop = FALSE
    ]
  )
  
  contrast_presence$group_role <- vapply(
    contrast_presence$group,
    group_role,
    character(1)
  )
  
  # Presenza nei due contrasti indipendentemente dal gruppo.
  contrast_term_summary <- aggregate(
    contrast ~ ontology + ID + Description,
    data = contrast_presence,
    FUN = function(x) paste(
      sort(unique(x)),
      collapse = " | "
    )
  )
  
  names(contrast_term_summary)[
    names(contrast_term_summary) == "contrast"
  ] <- "significant_contrasts"
  
  contrast_term_summary$number_of_contrasts <- vapply(
    strsplit(
      contrast_term_summary$significant_contrasts,
      " \\| "
    ),
    length,
    integer(1)
  )
  
  contrast_term_summary$interpretation <- ifelse(
    contrast_term_summary$number_of_contrasts == 2,
    "SIGNIFICANT_IN_BOTH_KNOCKDOWNS",
    "SIGNIFICANT_IN_ONE_KNOCKDOWN"
  )
  
  write_csv_safe(
    contrast_term_summary,
    file.path(TABLE_DIR, "term_comparison_between_knockdowns.csv")
  )
  
  # Confronto più stringente:
  # termini significativi nel gruppo COMMON_ALL_METHODS.
  robust_terms <- contrast_presence[
    contrast_presence$group == "COMMON_ALL_METHODS",
    ,
    drop = FALSE
  ]
  
  if (nrow(robust_terms) > 0) {
    
    robust_comparison <- aggregate(
      contrast ~ ontology + ID + Description,
      data = robust_terms,
      FUN = function(x) paste(
        sort(unique(x)),
        collapse = " | "
      )
    )
    
    names(robust_comparison)[
      names(robust_comparison) == "contrast"
    ] <- "significant_contrasts"
    
    robust_comparison$number_of_contrasts <- vapply(
      strsplit(
        robust_comparison$significant_contrasts,
        " \\| "
      ),
      length,
      integer(1)
    )
    
    robust_comparison$interpretation <- ifelse(
      robust_comparison$number_of_contrasts == 2,
      "ROBUST_SHARED_SIGNAL",
      "ROBUST_SIGNAL_IN_ONE_KNOCKDOWN"
    )
    
  } else {
    
    robust_comparison <- data.frame()
  }
  
  write_csv_safe(
    robust_comparison,
    file.path(TABLE_DIR, "robust_terms_between_knockdowns.csv")
  )
}


# ------------------------------------------------------------
# 10. CONFRONTO DELLE FIRME GENICHE DEI TERMINI
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("CONFRONTO DELLE FIRME GENICHE\n")
cat("------------------------------------------------------------\n\n")


if (nrow(significant_terms) > 0) {
  
  gene_signature_summary <- significant_terms[
    ,
    c(
      "contrast",
      "group",
      "ontology",
      "ID",
      "Description",
      "geneID"
    ),
    drop = FALSE
  ]
  
  gene_signature_summary$gene_count_parsed <- vapply(
    significant_terms$gene_list,
    length,
    integer(1)
  )
  
  gene_signature_summary$genes <- vapply(
    significant_terms$gene_list,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
  
  write_csv_safe(
    gene_signature_summary,
    file.path(TABLE_DIR, "significant_term_gene_signatures.csv")
  )
}


# ------------------------------------------------------------
# 11. RIEPILOGO DELLA CONSOLIDAZIONE
# ------------------------------------------------------------

if (nrow(consolidated_terms) > 0) {
  
  cluster_summary <- unique(
    consolidated_terms[
      ,
      c(
        "contrast",
        "group",
        "ontology",
        "cluster_id",
        "cluster_size",
        "ID",
        "Description",
        "p.adjust",
        "Count",
        "is_representative"
      ),
      drop = FALSE
    ]
  )
  
  cluster_summary <- cluster_summary[
    order(
      cluster_summary$contrast,
      cluster_summary$group,
      cluster_summary$ontology,
      cluster_summary$cluster_id,
      !cluster_summary$is_representative
    ),
    ,
    drop = FALSE
  ]
  
  write_csv_safe(
    cluster_summary,
    file.path(TABLE_DIR, "functional_clusters_summary.csv")
  )
  
  consolidation_summary <- aggregate(
    is_representative ~ contrast + group + ontology,
    data = consolidated_terms,
    FUN = function(x) {
      c(
        significant_terms = length(x),
        consolidated_clusters = sum(x)
      )
    }
  )
  
  consolidation_summary <- do.call(
    data.frame,
    consolidation_summary
  )
  
  names(consolidation_summary) <- c(
    "contrast",
    "group",
    "ontology",
    "significant_terms",
    "consolidated_clusters"
  )
  
  consolidation_summary$terms_reduced <- (
    consolidation_summary$significant_terms -
      consolidation_summary$consolidated_clusters
  )
  
  write_csv_safe(
    consolidation_summary,
    file.path(TABLE_DIR, "consolidation_summary.csv")
  )
  
} else {
  
  consolidation_summary <- data.frame()
}


# ------------------------------------------------------------
# 12. REPORT TESTUALE AUTOMATICO
# ------------------------------------------------------------

report_file <- file.path(
  REPORT_DIR,
  "functional_enrichment_interpretation_report.txt"
)

report_lines <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "FUNCTIONAL ENRICHMENT INTERPRETATION",
  "============================================================",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Input directory:", INPUT_DIR),
  paste("FDR threshold:", FDR_THRESHOLD),
  paste("Jaccard threshold:", JACCARD_THRESHOLD),
  paste("Containment threshold:", CONTAINMENT_THRESHOLD),
  "",
  "METHODOLOGICAL NOTES",
  "------------------------------------------------------------",
  "The analysis uses the enrichment CSV files generated by Script 31.",
  "No enrichment analysis has been recalculated.",
  "Significance is defined as BH-adjusted p-value < 0.05.",
  "Term consolidation is based on overlap between associated gene IDs.",
  "It does not represent semantic GO similarity.",
  "The representative term is selected by lowest adjusted p-value,",
  "with Count used as a secondary ordering criterion.",
  "",
  "IMPORT SUMMARY",
  "------------------------------------------------------------",
  paste("CSV files discovered:", length(csv_files)),
  paste("Valid enrichment files:", length(imported_results)),
  paste("Imported enrichment records:", nrow(enrichment)),
  paste("Significant terms:", nrow(significant_terms)),
  "",
  "SIGNIFICANT TERMS BY CONTRAST, GROUP AND ONTOLOGY",
  "------------------------------------------------------------"
)

if (nrow(enrichment_summary) > 0) {
  
  for (i in seq_len(nrow(enrichment_summary))) {
    
    current <- enrichment_summary[i, ]
    
    report_lines <- c(
      report_lines,
      paste0(
        current$contrast,
        " | ",
        current$group,
        " | ",
        current$ontology,
        " | total=",
        current$total_terms,
        " | significant=",
        current$significant_terms
      )
    )
  }
}

report_lines <- c(
  report_lines,
  "",
  "CONSOLIDATION",
  "------------------------------------------------------------"
)

if (nrow(consolidation_summary) > 0) {
  
  for (i in seq_len(nrow(consolidation_summary))) {
    
    current <- consolidation_summary[i, ]
    
    report_lines <- c(
      report_lines,
      paste0(
        current$contrast,
        " | ",
        current$group,
        " | ",
        current$ontology,
        " | significant terms=",
        current$significant_terms,
        " | clusters=",
        current$consolidated_clusters,
        " | terms reduced=",
        current$terms_reduced
      )
    )
  }
  
} else {
  
  report_lines <- c(
    report_lines,
    "No significant terms available for consolidation."
  )
}

report_lines <- c(
  report_lines,
  "",
  "INTERPRETATION GUIDANCE",
  "------------------------------------------------------------",
  "COMMON_ALL_METHODS:",
  "Terms significant across all normalization methods; strongest",
  "candidate signals for normalization robustness.",
  "",
  "DISCORDANT_ANY_METHOD:",
  "Terms significant in at least one method but not consistently",
  "supported by all methods.",
  "",
  "ALTERNATIVE_ONLY:",
  "Terms not significant with TMM but significant with at least",
  "one alternative normalization method.",
  "",
  "MIXED_METHOD_DISCORDANCE:",
  "Terms whose significance varies among normalization methods.",
  "",
  "EXCLUSIVE_*:",
  "Terms associated with genes significant only under a specific",
  "alternative normalization method.",
  "",
  "IMPORTANT:",
  "Enrichment significance does not demonstrate pathway activation,",
  "causality, or a direct functional effect of NRP1 knockdown.",
  "Interpretation requires examination of term descriptions,",
  "gene membership, direction of expression and biological context.",
  "",
  "============================================================",
  "END OF REPORT",
  "============================================================"
)

writeLines(
  report_lines,
  con = report_file,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 13. AUDIT FINALE
# ------------------------------------------------------------

audit_summary <- data.frame(
  metric = c(
    "csv_files_discovered",
    "valid_enrichment_files",
    "imported_enrichment_records",
    "significant_terms",
    "consolidated_terms",
    "representative_terms",
    "contrasts",
    "groups",
    "ontologies"
  ),
  value = c(
    length(csv_files),
    length(imported_results),
    nrow(enrichment),
    nrow(significant_terms),
    nrow(consolidated_terms),
    nrow(representative_terms),
    length(unique(enrichment$contrast)),
    length(unique(enrichment$group)),
    length(unique(enrichment$ontology))
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  audit_summary,
  file.path(AUDIT_DIR, "final_audit_summary.csv")
)

# Salvataggio degli oggetti analitici per eventuali analisi successive.
saveRDS(
  list(
    enrichment = enrichment,
    significant_terms = significant_terms,
    consolidated_terms = consolidated_terms,
    representative_terms = representative_terms,
    enrichment_summary = enrichment_summary,
    consolidation_summary = consolidation_summary,
    import_log = import_log_df,
    parameters = list(
      FDR_THRESHOLD = FDR_THRESHOLD,
      JACCARD_THRESHOLD = JACCARD_THRESHOLD,
      CONTAINMENT_THRESHOLD = CONTAINMENT_THRESHOLD
    )
  ),
  file = file.path(
    OUTPUT_DIR,
    "functional_enrichment_interpretation_objects.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    AUDIT_DIR,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 14. OUTPUT FINALE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FUNCTIONAL ENRICHMENT INTERPRETATION COMPLETATA\n")
cat("============================================================\n\n")

cat("Record di arricchimento:", nrow(enrichment), "\n")
cat("Termini significativi:", nrow(significant_terms), "\n")

if (nrow(consolidated_terms) > 0) {
  
  cat(
    "Cluster funzionali:",
    length(unique(consolidated_terms$cluster_id)),
    "\n"
  )
  
  cat(
    "Termini rappresentativi:",
    nrow(representative_terms),
    "\n"
  )
  
} else {
  
  cat("Cluster funzionali: 0\n")
  cat("Termini rappresentativi: 0\n")
}

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n\n")

cat("File principali:\n")
cat("- complete_enrichment_inventory.csv\n")
cat("- significant_enrichment_terms.csv\n")
cat("- redundancy_consolidated_terms.csv\n")
cat("- representative_functional_terms.csv\n")
cat("- functional_clusters_summary.csv\n")
cat("- term_comparison_between_knockdowns.csv\n")
cat("- robust_terms_between_knockdowns.csv\n")
cat("- significant_term_gene_signatures.csv\n")
cat("- functional_enrichment_interpretation_report.txt\n\n")

cat("Nessun risultato originale è stato modificato.\n")
cat("============================================================\n")
