
# ============================================================
# TRANSCRIPTOME-X
# Script: 31_discordant_gene_functional_analysis.R
#
# DISCORDANT GENE FUNCTIONAL ANALYSIS
#
# Versione robusta con:
# - gestione dei gruppi vuoti;
# - recupero dei risultati CSV gia' prodotti;
# - cache per singola analisi;
# - mapping Ensembl -> Entrez;
# - enrichment GO BP, CC, MF e Reactome;
# - universo composto da tutti i geni testati e mappabili;
# - gestione degli errori per singola analisi;
# - esportazione CSV, PDF, RDS e report;
# - controlli di coerenza.
#
# Non ricalcola normalizzazione o modelli DE.
# ============================================================


# ------------------------------------------------------------
# 1. CONFIGURAZIONE
# ------------------------------------------------------------

project_dir <- "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X"

input_dir <- file.path(
  project_dir,
  "results",
  "normalization_discordant_genes"
)

output_dir <- file.path(
  project_dir,
  "results",
  "discordant_gene_functional_analysis"
)

contrasts <- c(
  "siNRP1_3_vs_siNT",
  "siNRP1_5_vs_siNT"
)

fdr_threshold <- 0.05
min_gene_set_size <- 10
max_gene_set_size <- 500
top_terms_plot <- 15

ontologies <- c(
  "BP",
  "CC",
  "MF",
  "Reactome"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (contrast in contrasts) {
  dir.create(
    file.path(output_dir, contrast),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 2. DIPENDENZE
# ------------------------------------------------------------

required_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ReactomePA"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Pacchetti mancanti: ",
    paste(missing_packages, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 3. INTESTAZIONE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("DISCORDANT GENE FUNCTIONAL ANALYSIS\n")
cat("============================================================\n\n")

cat("Pacchetti verificati.\n")

cat(
  "clusterProfiler:",
  as.character(utils::packageVersion("clusterProfiler")),
  "\n"
)

cat(
  "ReactomePA:",
  as.character(utils::packageVersion("ReactomePA")),
  "\n"
)


# ------------------------------------------------------------
# 4. VERIFICA INPUT
# ------------------------------------------------------------

if (!dir.exists(input_dir)) {
  stop("Directory input non trovata: ", input_dir)
}

input_files <- file.path(
  input_dir,
  paste0(
    "discordant_gene_audit_",
    contrasts,
    ".csv"
  )
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "File di input mancanti:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("\nDirectory input verificata.\n")
cat("Contrasti disponibili:", length(contrasts), "\n")


# ------------------------------------------------------------
# 5. FUNZIONI GENERALI
# ------------------------------------------------------------

write_csv <- function(data, filename) {
  
  if (!is.data.frame(data)) {
    stop("write_csv richiede un data.frame.")
  }
  
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  temporary_file <- tempfile(
    pattern = "transcriptome_x_",
    tmpdir = dirname(filename),
    fileext = ".csv"
  )
  
  utils::write.csv(
    data,
    file = temporary_file,
    row.names = FALSE,
    na = ""
  )
  
  if (!file.rename(temporary_file, filename)) {
    if (file.exists(temporary_file)) {
      unlink(temporary_file)
    }
    
    stop("Impossibile salvare il file: ", filename)
  }
  
  invisible(filename)
}


empty_enrichment_table <- function() {
  
  data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
}


read_audit <- function(contrast) {
  
  filepath <- file.path(
    input_dir,
    paste0(
      "discordant_gene_audit_",
      contrast,
      ".csv"
    )
  )
  
  data <- utils::read.csv(
    filepath,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "discordance_class",
    "significant_methods"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Colonne mancanti in ",
      basename(filepath),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data$gene_id <- trimws(
    as.character(data$gene_id)
  )
  
  data$discordance_class <- as.character(
    data$discordance_class
  )
  
  if (
    anyNA(data$gene_id) ||
    any(data$gene_id == "")
  ) {
    stop("Gene ID mancanti nel file: ", basename(filepath))
  }
  
  if (anyDuplicated(data$gene_id) > 0) {
    stop("Gene ID duplicati nel file: ", basename(filepath))
  }
  
  if (anyNA(data$discordance_class)) {
    stop("Classi di discordanza mancanti: ", basename(filepath))
  }
  
  data
}


# ------------------------------------------------------------
# 6. MAPPING ENSEMBL -> ENTREZ
# ------------------------------------------------------------

map_gene_ids <- function(gene_ids) {
  
  original_ids <- unique(
    trimws(as.character(gene_ids))
  )
  
  clean_ids <- sub(
    "\\.[0-9]+$",
    "",
    original_ids
  )
  
  clean_ids[clean_ids == ""] <- NA_character_
  
  valid_ids <- unique(
    clean_ids[!is.na(clean_ids)]
  )
  
  if (length(valid_ids) == 0) {
    stop("Nessun Ensembl ID valido da mappare.")
  }
  
  entrez_mapping <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_ids,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  symbol_mapping <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  lookup <- data.frame(
    Ensembl_ID = valid_ids,
    ENTREZID = unname(
      as.character(entrez_mapping[valid_ids])
    ),
    gene_symbol = unname(
      as.character(symbol_mapping[valid_ids])
    ),
    stringsAsFactors = FALSE
  )
  
  index <- match(
    clean_ids,
    lookup$Ensembl_ID
  )
  
  result <- data.frame(
    gene_id = original_ids,
    Ensembl_ID = clean_ids,
    ENTREZID = lookup$ENTREZID[index],
    gene_symbol = lookup$gene_symbol[index],
    stringsAsFactors = FALSE
  )
  
  result$mapped_entrez <- (
    !is.na(result$ENTREZID) &
      result$ENTREZID != ""
  )
  
  result$mapped_symbol <- (
    !is.na(result$gene_symbol) &
      result$gene_symbol != ""
  )
  
  result
}


# ------------------------------------------------------------
# 7. DEFINIZIONE DEI GRUPPI
# ------------------------------------------------------------

build_gene_groups <- function(data) {
  
  groups <- list()
  
  get_class <- function(class_name) {
    
    data$gene_id[
      data$discordance_class == class_name
    ]
  }
  
  get_exclusive <- function(column_name) {
    
    if (!column_name %in% colnames(data)) {
      return(character())
    }
    
    values <- data[[column_name]]
    
    if (is.logical(values)) {
      
      selected <- !is.na(values) & values
      
    } else {
      
      values <- tolower(trimws(as.character(values)))
      
      selected <- !is.na(values) &
        values %in% c("true", "t", "1", "yes")
    }
    
    data$gene_id[selected]
  }
  
  groups[["COMMON_ALL_METHODS"]] <- get_class(
    "SIGNIFICANT_ALL_METHODS"
  )
  
  groups[["DISCORDANT_ANY_METHOD"]] <- data$gene_id[
    data$discordance_class %in% c(
      "TMM_ONLY",
      "ALTERNATIVE_ONLY",
      "MIXED_METHOD_DISCORDANCE"
    )
  ]
  
  groups[["ALTERNATIVE_ONLY"]] <- get_class(
    "ALTERNATIVE_ONLY"
  )
  
  groups[["MIXED_METHOD_DISCORDANCE"]] <- get_class(
    "MIXED_METHOD_DISCORDANCE"
  )
  
  groups[["TMM_ONLY"]] <- get_class(
    "TMM_ONLY"
  )
  
  groups[["EXCLUSIVE_TMMwsp"]] <- get_exclusive(
    "significant_only_TMMwsp"
  )
  
  groups[["EXCLUSIVE_RLE"]] <- get_exclusive(
    "significant_only_RLE"
  )
  
  groups[["EXCLUSIVE_UPPERQUARTILE"]] <- get_exclusive(
    "significant_only_upperquartile"
  )
  
  groups <- lapply(
    groups,
    function(x) {
      unique(as.character(x))
    }
  )
  
  groups
}


# ------------------------------------------------------------
# 8. ENRICHMENT
# ------------------------------------------------------------

run_enrichment <- function(
    gene_entrez,
    universe_entrez,
    ontology
) {
  
  gene_entrez <- unique(
    as.character(gene_entrez)
  )
  
  universe_entrez <- unique(
    as.character(universe_entrez)
  )
  
  gene_entrez <- gene_entrez[
    !is.na(gene_entrez) &
      gene_entrez != ""
  ]
  
  universe_entrez <- universe_entrez[
    !is.na(universe_entrez) &
      universe_entrez != ""
  ]
  
  gene_entrez <- intersect(
    gene_entrez,
    universe_entrez
  )
  
  if (length(gene_entrez) == 0) {
    return(list(
      data = empty_enrichment_table(),
      status = "SKIPPED_EMPTY_GROUP",
      message = "Nessun gene Entrez nel gruppo."
    ))
  }
  
  if (length(gene_entrez) < min_gene_set_size) {
    return(list(
      data = empty_enrichment_table(),
      status = "SKIPPED_TOO_FEW_GENES",
      message = paste(
        "Geni mappati:",
        length(gene_entrez),
        "- minimo richiesto:",
        min_gene_set_size
      )
    ))
  }
  
  if (length(universe_entrez) < min_gene_set_size) {
    return(list(
      data = empty_enrichment_table(),
      status = "SKIPPED_UNIVERSE_TOO_SMALL",
      message = "Universo di riferimento insufficiente."
    ))
  }
  
  enrichment_result <- tryCatch(
    
    {
      
      if (ontology %in% c("BP", "CC", "MF")) {
        
        clusterProfiler::enrichGO(
          gene = gene_entrez,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = ontology,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = min_gene_set_size,
          maxGSSize = max_gene_set_size,
          readable = TRUE
        )
        
      } else if (ontology == "Reactome") {
        
        ReactomePA::enrichPathway(
          gene = gene_entrez,
          universe = universe_entrez,
          organism = "human",
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = min_gene_set_size,
          maxGSSize = max_gene_set_size,
          readable = TRUE
        )
        
      } else {
        
        stop("Ontologia non riconosciuta: ", ontology)
      }
    },
    
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "enrichment_error"
      )
    }
  )
  
  if (inherits(enrichment_result, "enrichment_error")) {
    return(list(
      data = empty_enrichment_table(),
      status = "ERROR",
      message = enrichment_result$message
    ))
  }
  
  result_table <- tryCatch(
    as.data.frame(enrichment_result),
    error = function(e) NULL
  )
  
  if (
    is.null(result_table) ||
    !is.data.frame(result_table) ||
    nrow(result_table) == 0
  ) {
    return(list(
      data = empty_enrichment_table(),
      status = "NO_ENRICHMENT",
      message = "Nessun termine restituito."
    ))
  }
  
  if (!"p.adjust" %in% colnames(result_table)) {
    return(list(
      data = empty_enrichment_table(),
      status = "ERROR",
      message = "Colonna p.adjust assente nei risultati."
    ))
  }
  
  result_table <- as.data.frame(
    result_table,
    stringsAsFactors = FALSE
  )
  
  list(
    data = result_table,
    status = "COMPLETED",
    message = paste(
      nrow(result_table),
      "termini restituiti."
    )
  )
}


# ------------------------------------------------------------
# 9. RECUPERO RISULTATI ESISTENTI
# ------------------------------------------------------------

load_existing_result <- function(filepath) {
  
  if (!file.exists(filepath)) {
    return(NULL)
  }
  
  if (file.info(filepath)$size <= 0) {
    return(NULL)
  }
  
  cached <- tryCatch(
    utils::read.csv(
      filepath,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(cached) ||
    !is.data.frame(cached) ||
    nrow(cached) == 0 ||
    !"p.adjust" %in% colnames(cached) ||
    !"Description" %in% colnames(cached) ||
    !"ID" %in% colnames(cached)
  ) {
    return(NULL)
  }
  
  cached
}


# ------------------------------------------------------------
# 10. GRAFICI
# ------------------------------------------------------------

plot_enrichment_pdf <- function(
    results_list,
    contrast,
    output_file
) {
  
  grDevices::pdf(
    output_file,
    width = 12,
    height = 9,
    onefile = TRUE
  )
  
  on.exit(
    grDevices::dev.off(),
    add = TRUE
  )
  
  for (ontology in ontologies) {
    
    combined <- list()
    
    for (group_name in names(results_list)) {
      
      current <- results_list[[group_name]][[ontology]]$data
      
      if (
        is.null(current) ||
        !is.data.frame(current) ||
        nrow(current) == 0
      ) {
        next
      }
      
      if (!"p.adjust" %in% colnames(current)) {
        next
      }
      
      current <- current[
        !is.na(current$p.adjust) &
          is.finite(current$p.adjust),
        ,
        drop = FALSE
      ]
      
      if (nrow(current) == 0) {
        next
      }
      
      current$group <- group_name
      
      combined[[group_name]] <- current
    }
    
    if (length(combined) == 0) {
      
      graphics::plot.new()
      
      graphics::title(
        main = paste(contrast, "-", ontology)
      )
      
      graphics::text(
        0.5,
        0.5,
        "Nessun arricchimento disponibile."
      )
      
      next
    }
    
    plot_data <- do.call(
      rbind,
      combined
    )
    
    plot_data <- plot_data[
      order(plot_data$p.adjust),
      ,
      drop = FALSE
    ]
    
    plot_data <- utils::head(
      plot_data,
      top_terms_plot
    )
    
    plot_data$plot_label <- paste(
      plot_data$group,
      "-",
      plot_data$Description
    )
    
    plot_data$minus_log10_FDR <- -log10(
      pmax(
        plot_data$p.adjust,
        .Machine$double.xmin
      )
    )
    
    graphics::par(
      mar = c(5, 18, 4, 2)
    )
    
    graphics::barplot(
      rev(plot_data$minus_log10_FDR),
      names.arg = rev(plot_data$plot_label),
      horiz = TRUE,
      las = 1,
      cex.names = 0.65,
      xlab = "-log10 adjusted p-value",
      main = paste(
        contrast,
        "-",
        ontology,
        "(top terms)"
      )
    )
    
    graphics::abline(
      v = -log10(fdr_threshold),
      lty = 2
    )
  }
  
  invisible(NULL)
}


# ------------------------------------------------------------
# 11. FUNZIONE DI MAPPING DEI GRUPPI
# ------------------------------------------------------------

map_group_to_entrez <- function(gene_ids, mapping) {
  
  if (length(gene_ids) == 0) {
    return(character())
  }
  
  selected <- mapping[
    mapping$gene_id %in% gene_ids &
      mapping$mapped_entrez,
    ,
    drop = FALSE
  ]
  
  result <- unique(
    as.character(selected$ENTREZID)
  )
  
  result <- result[
    !is.na(result) &
      result != ""
  ]
  
  result
}


# ------------------------------------------------------------
# 12. ANALISI PER CONTRASTO
# ------------------------------------------------------------

all_status <- list()
all_mapping_audits <- list()
all_group_summaries <- list()
all_enrichment_results <- list()

for (contrast in contrasts) {
  
  cat("\n")
  cat("============================================================\n")
  cat("CONTRASTO:", contrast, "\n")
  cat("============================================================\n")
  
  contrast_output <- file.path(
    output_dir,
    contrast
  )
  
  audit <- read_audit(contrast)
  
  cat(
    "Geni nell'audit:",
    nrow(audit),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 12A. MAPPING
  # ----------------------------------------------------------
  
  mapping <- map_gene_ids(
    audit$gene_id
  )
  
  if (nrow(mapping) != nrow(audit)) {
    stop("Incoerenza nel numero di righe del mapping.")
  }
  
  if (!identical(mapping$gene_id, audit$gene_id)) {
    stop("L'ordine dei gene ID nel mapping non corrisponde all'audit.")
  }
  
  audit_mapped <- audit
  
  audit_mapped$Ensembl_ID_mapping <- mapping$Ensembl_ID
  audit_mapped$ENTREZID <- mapping$ENTREZID
  audit_mapped$gene_symbol_mapping <- mapping$gene_symbol
  audit_mapped$mapped_entrez <- mapping$mapped_entrez
  audit_mapped$mapped_symbol <- mapping$mapped_symbol
  
  mapping_summary <- data.frame(
    contrast = contrast,
    genes_tested = nrow(audit),
    genes_with_entrez = sum(mapping$mapped_entrez),
    genes_without_entrez = sum(!mapping$mapped_entrez),
    genes_with_symbol = sum(mapping$mapped_symbol),
    unique_entrez_ids = length(
      unique(
        mapping$ENTREZID[mapping$mapped_entrez]
      )
    ),
    stringsAsFactors = FALSE
  )
  
  all_mapping_audits[[contrast]] <- mapping_summary
  
  write_csv(
    audit_mapped,
    file.path(
      contrast_output,
      "gene_mapping_audit.csv"
    )
  )
  
  universe_entrez <- unique(
    mapping$ENTREZID[mapping$mapped_entrez]
  )
  
  universe_entrez <- universe_entrez[
    !is.na(universe_entrez) &
      universe_entrez != ""
  ]
  
  if (length(universe_entrez) == 0) {
    stop("Nessun gene mappato nell'universo per ", contrast)
  }
  
  cat(
    "Geni con Entrez ID:",
    mapping_summary$genes_with_entrez,
    "\n"
  )
  
  cat(
    "Universo Entrez unici:",
    length(universe_entrez),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 12B. GRUPPI
  # ----------------------------------------------------------
  
  gene_groups <- build_gene_groups(audit)
  
  group_mapping <- lapply(
    gene_groups,
    map_group_to_entrez,
    mapping = mapping
  )
  
  group_summary <- data.frame(
    contrast = contrast,
    group = names(gene_groups),
    genes_original = vapply(
      gene_groups,
      length,
      integer(1)
    ),
    genes_mapped = vapply(
      group_mapping,
      length,
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  
  group_summary$percentage_mapped <- ifelse(
    group_summary$genes_original > 0,
    100 * group_summary$genes_mapped /
      group_summary$genes_original,
    NA_real_
  )
  
  all_group_summaries[[contrast]] <- group_summary
  
  write_csv(
    group_summary,
    file.path(
      contrast_output,
      "gene_group_mapping_summary.csv"
    )
  )
  
  for (group_name in names(gene_groups)) {
    
    current_ids <- gene_groups[[group_name]]
    
    current_data <- audit_mapped[
      audit_mapped$gene_id %in% current_ids,
      ,
      drop = FALSE
    ]
    
    write_csv(
      current_data,
      file.path(
        contrast_output,
        paste0(
          "gene_list_",
          group_name,
          ".csv"
        )
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 12C. ENRICHMENT GO E REACTOME
  # ----------------------------------------------------------
  
  contrast_results <- list()
  contrast_status <- list()
  
  for (group_name in names(gene_groups)) {
    
    cat("\n")
    cat("Gruppo:", group_name, "\n")
    
    contrast_results[[group_name]] <- list()
    
    for (ontology in ontologies) {
      
      output_name <- paste0(
        "enrichment_",
        group_name,
        "_",
        ontology,
        ".csv"
      )
      
      output_path <- file.path(
        contrast_output,
        output_name
      )
      
      # ------------------------------------------------------
      # Recupero di un CSV gia' prodotto
      # ------------------------------------------------------
      
      cached_result <- load_existing_result(
        output_path
      )
      
      if (!is.null(cached_result)) {
        
        result <- list(
          data = cached_result,
          status = "RECOVERED_EXISTING_CSV",
          message = paste(
            "Risultato recuperato dal CSV esistente:",
            nrow(cached_result),
            "termini."
          )
        )
        
        cat(
          ontology,
          ": RECOVERED_EXISTING_CSV -",
          nrow(cached_result),
          "termini\n"
        )
        
      } else {
        
        # ----------------------------------------------------
        # Nuova analisi
        # ----------------------------------------------------
        
        result <- run_enrichment(
          gene_entrez = group_mapping[[group_name]],
          universe_entrez = universe_entrez,
          ontology = ontology
        )
        
        n_results <- nrow(result$data)
        
        result$data$contrast <- rep(
          contrast,
          n_results
        )
        
        result$data$group <- rep(
          group_name,
          n_results
        )
        
        result$data$ontology <- rep(
          ontology,
          n_results
        )
        
        write_csv(
          result$data,
          output_path
        )
        
        cat(
          ontology,
          ":",
          result$status,
          "-",
          nrow(result$data),
          "termini\n"
        )
      }
      
      # I risultati recuperati hanno gia' le colonne
      # contrast, group e ontology.
      if (nrow(result$data) > 0) {
        
        if (!"contrast" %in% colnames(result$data)) {
          result$data$contrast <- rep(
            contrast,
            nrow(result$data)
          )
        }
        
        if (!"group" %in% colnames(result$data)) {
          result$data$group <- rep(
            group_name,
            nrow(result$data)
          )
        }
        
        if (!"ontology" %in% colnames(result$data)) {
          result$data$ontology <- rep(
            ontology,
            nrow(result$data)
          )
        }
      }
      
      contrast_results[[group_name]][[ontology]] <- result
      
      significant_count <- if (
        nrow(result$data) > 0 &&
        "p.adjust" %in% colnames(result$data)
      ) {
        sum(
          !is.na(result$data$p.adjust) &
            result$data$p.adjust < fdr_threshold
        )
      } else {
        0L
      }
      
      status_row <- data.frame(
        contrast = contrast,
        group = group_name,
        ontology = ontology,
        genes_original = length(
          gene_groups[[group_name]]
        ),
        genes_mapped = length(
          group_mapping[[group_name]]
        ),
        universe_size = length(universe_entrez),
        status = result$status,
        terms_returned = nrow(result$data),
        significant_terms = significant_count,
        message = result$message,
        stringsAsFactors = FALSE
      )
      
      contrast_status[[length(contrast_status) + 1L]] <-
        status_row
    }
  }
  
  contrast_status_table <- do.call(
    rbind,
    contrast_status
  )
  
  rownames(contrast_status_table) <- NULL
  
  all_status[[contrast]] <- contrast_status_table
  
  write_csv(
    contrast_status_table,
    file.path(
      contrast_output,
      "enrichment_execution_summary.csv"
    )
  )
  
  
  # ----------------------------------------------------------
  # 12D. RISULTATI SIGNIFICATIVI
  # ----------------------------------------------------------
  
  significant_results <- list()
  
  for (group_name in names(contrast_results)) {
    
    for (ontology in ontologies) {
      
      current <- contrast_results[[group_name]][[ontology]]$data
      
      if (
        is.null(current) ||
        !is.data.frame(current) ||
        nrow(current) == 0 ||
        !"p.adjust" %in% colnames(current)
      ) {
        next
      }
      
      current <- current[
        !is.na(current$p.adjust) &
          current$p.adjust < fdr_threshold,
        ,
        drop = FALSE
      ]
      
      if (nrow(current) > 0) {
        
        significant_results[[
          paste(group_name, ontology, sep = "__")
        ]] <- current
      }
    }
  }
  
  if (length(significant_results) > 0) {
    
    significant_table <- do.call(
      rbind,
      significant_results
    )
    
    rownames(significant_table) <- NULL
    
  } else {
    
    significant_table <- data.frame(
      contrast = character(),
      group = character(),
      ontology = character(),
      ID = character(),
      Description = character(),
      p.adjust = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  write_csv(
    significant_table,
    file.path(
      contrast_output,
      "all_significant_enrichment_terms.csv"
    )
  )
  
  
  # ----------------------------------------------------------
  # 12E. GRAFICI
  # ----------------------------------------------------------
  
  plot_enrichment_pdf(
    results_list = contrast_results,
    contrast = contrast,
    output_file = file.path(
      contrast_output,
      "functional_enrichment_top_terms.pdf"
    )
  )
  
  
  # ----------------------------------------------------------
  # 12F. OGGETTI RDS
  # ----------------------------------------------------------
  
  saveRDS(
    list(
      contrast = contrast,
      audit = audit,
      mapping = mapping,
      universe_entrez = universe_entrez,
      gene_groups = gene_groups,
      group_mapping = group_mapping,
      group_summary = group_summary,
      enrichment_results = contrast_results,
      execution_summary = contrast_status_table,
      parameters = list(
        fdr_threshold = fdr_threshold,
        min_gene_set_size = min_gene_set_size,
        max_gene_set_size = max_gene_set_size,
        ontologies = ontologies
      )
    ),
    file = file.path(
      contrast_output,
      "functional_analysis_complete.rds"
    )
  )
  
  all_enrichment_results[[contrast]] <- contrast_results
  
  cat("\n")
  cat("Contrasto completato:", contrast, "\n")
}


# ------------------------------------------------------------
# 13. RIEPILOGHI GLOBALI
# ------------------------------------------------------------

mapping_summary_all <- do.call(
  rbind,
  all_mapping_audits
)

group_summary_all <- do.call(
  rbind,
  all_group_summaries
)

status_summary_all <- do.call(
  rbind,
  all_status
)

rownames(mapping_summary_all) <- NULL
rownames(group_summary_all) <- NULL
rownames(status_summary_all) <- NULL

write_csv(
  mapping_summary_all,
  file.path(
    output_dir,
    "global_mapping_summary.csv"
  )
)

write_csv(
  group_summary_all,
  file.path(
    output_dir,
    "global_gene_group_summary.csv"
  )
)

write_csv(
  status_summary_all,
  file.path(
    output_dir,
    "global_enrichment_execution_summary.csv"
  )
)


# ------------------------------------------------------------
# 14. REPORT FINALE
# ------------------------------------------------------------

report_path <- file.path(
  output_dir,
  "DISCORDANT_GENE_FUNCTIONAL_ANALYSIS_REPORT.txt"
)

report_lines <- c(
  "============================================================",
  "TRANSCRIPTOME-X",
  "DISCORDANT GENE FUNCTIONAL ANALYSIS",
  "FINAL REPORT",
  "============================================================",
  "",
  "PARAMETRI ANALISI",
  "------------------------------------------------------------",
  paste("Soglia FDR:", fdr_threshold),
  paste("Min gene set size:", min_gene_set_size),
  paste("Max gene set size:", max_gene_set_size),
  "Correzione multipli test: Benjamini-Hochberg",
  "Organismo: Homo sapiens",
  "Universo: tutti i geni testati e mappabili a Entrez ID",
  "Database: org.Hs.eg.db",
  "Ontologie: GO BP, GO CC, GO MF, Reactome",
  "",
  "MAPPING SUMMARY",
  "------------------------------------------------------------",
  capture.output(
    print(mapping_summary_all, row.names = FALSE)
  ),
  "",
  "GENE GROUP SUMMARY",
  "------------------------------------------------------------",
  capture.output(
    print(group_summary_all, row.names = FALSE)
  ),
  "",
  "ENRICHMENT EXECUTION SUMMARY",
  "------------------------------------------------------------",
  capture.output(
    print(status_summary_all, row.names = FALSE)
  ),
  "",
  "INTERPRETAZIONE METODOLOGICA",
  "------------------------------------------------------------",
  paste(
    "L'enrichment e' stato eseguito separatamente",
    "per ciascun contrasto e gruppo."
  ),
  paste(
    "L'universo e' costituito dai geni testati",
    "e mappabili a Entrez ID."
  ),
  paste(
    "Le liste con meno di",
    min_gene_set_size,
    "geni mappati non sono state analizzate."
  ),
  paste(
    "I risultati sono stati corretti con Benjamini-Hochberg",
    "all'interno di ciascuna analisi."
  ),
  "Un risultato non significativo non dimostra assenza di rilevanza biologica.",
  "L'enrichment non dimostra attivazione o inibizione di pathway.",
  "Le liste esclusive dei metodi possono avere potenza limitata.",
  "I termini arricchiti non dimostrano effetti causali del knockdown.",
  paste(
    "L'analisi non risolve la limitazione relativa agli",
    "effective lengths RSEM non disponibili."
  ),
  "",
  "NOTA SUI RISULTATI RECUPERATI",
  "------------------------------------------------------------",
  paste(
    "I CSV validi gia' presenti sono stati recuperati",
    "senza ripetere le rispettive analisi."
  ),
  paste(
    "I gruppi vuoti sono stati registrati come",
    "SKIPPED_EMPTY_GROUP."
  )
)

writeLines(
  report_lines,
  con = report_path,
  useBytes = TRUE
)


# ------------------------------------------------------------
# 15. SESSION INFO
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 16. INVENTARIO FINALE
# ------------------------------------------------------------

generated_files <- list.files(
  output_dir,
  recursive = TRUE,
  full.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("DISCORDANT GENE FUNCTIONAL ANALYSIS - FINAL REPORT\n")
cat("============================================================\n\n")

cat("Mapping summary:\n")
print(
  mapping_summary_all,
  row.names = FALSE
)

cat("\nGene group summary:\n")
print(
  group_summary_all,
  row.names = FALSE
)

cat("\nEnrichment execution summary:\n")
print(
  status_summary_all,
  row.names = FALSE
)

cat("\nDirectory output:\n")
cat(
  normalizePath(output_dir),
  "\n"
)

cat(
  "\nNumero file generati:",
  length(generated_files),
  "\n"
)

cat("\nFile generati:\n")
print(generated_files)

cat("\n")
cat("ANALISI FUNZIONALE COMPLETATA.\n")
8