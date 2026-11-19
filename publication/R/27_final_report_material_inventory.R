
# ============================================================
# TRANSCRIPTOME-X
# Script: 27_final_report_material_inventory.R
#
# Obiettivo:
# Recuperare i risultati GO/Reactome e le figure gia' prodotti.
# Non esegue nuove analisi e non modifica gli output esistenti.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X - FINAL REPORT MATERIAL INVENTORY\n")
cat("============================================================\n\n")

project_dir <- getwd()

if (!dir.exists(file.path(project_dir, "data", "processed",
                          "functional_enrichment_strict"))) {
  stop("Cartella functional_enrichment_strict non trovata. ",
       "Controlla la working directory.")
}

base_dir <- file.path(
  project_dir,
  "data", "processed", "functional_enrichment_strict"
)

results_dir <- file.path(base_dir, "results")
plots_dir <- file.path(base_dir, "plots")
reports_dir <- file.path(base_dir, "reports")

output_dir <- file.path(
  reports_dir,
  "final_report_material_inventory"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. INVENTARIO DEI FILE DI ENRICHMENT
# ------------------------------------------------------------

cat("1. FILE DI ENRICHMENT\n")
cat("------------------------------------------------------------\n")

if (!dir.exists(results_dir)) {
  stop("Cartella results non trovata.")
}

result_files <- list.files(
  results_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(result_files) == 0) {
  cat("Nessun CSV trovato nella cartella results.\n")
} else {
  
  result_inventory <- data.frame(
    file = basename(result_files),
    relative_path = substring(
      result_files,
      nchar(results_dir) + 2
    ),
    size_kb = round(
      file.info(result_files)$size / 1024,
      2
    ),
    stringsAsFactors = FALSE
  )
  
  print(result_inventory, row.names = FALSE)
  
  write.csv(
    result_inventory,
    file.path(output_dir, "R27_enrichment_file_inventory.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 2. LETTURA RISULTATI ESTRATTI
# ------------------------------------------------------------

cat("\n2. RISULTATI GO E REACTOME\n")
cat("------------------------------------------------------------\n")

pathway_results <- list()

for (f in result_files) {
  
  filename <- basename(f)
  
  if (!grepl(
    "BP|MF|CC|REACTOME",
    filename,
    ignore.case = TRUE
  )) {
    next
  }
  
  dat <- tryCatch(
    read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  
  cat("\nFILE:", filename, "\n")
  cat("Righe:", nrow(dat), "\n")
  cat("Colonne:", paste(names(dat), collapse = " | "), "\n")
  
  # Individua la colonna statistica disponibile
  p_candidates <- c(
    "p.adjust",
    "FDR",
    "padj",
    "qvalue",
    "pvalue",
    "PValue"
  )
  
  p_col <- p_candidates[
    p_candidates %in% names(dat)
  ]
  
  # Individua la colonna descrittiva
  description_candidates <- c(
    "Description",
    "description",
    "term",
    "Term",
    "pathway",
    "Pathway"
  )
  
  description_col <- description_candidates[
    description_candidates %in% names(dat)
  ]
  
  if (length(description_col) == 0) {
    description_col <- NA_character_
  } else {
    description_col <- description_col[1]
  }
  
  if (length(p_col) > 0) {
    
    p_col <- p_col[1]
    
    dat[[p_col]] <- suppressWarnings(
      as.numeric(dat[[p_col]])
    )
    
    dat <- dat[
      !is.na(dat[[p_col]]),
      ,
      drop = FALSE
    ]
    
    dat <- dat[order(dat[[p_col]]), , drop = FALSE]
    
    selected <- head(dat, 10)
    
    cat("\nTop risultati ordinati per", p_col, ":\n")
    
    if (!is.na(description_col)) {
      
      display_cols <- unique(c(
        description_col,
        p_col,
        intersect(
          c("ID", "Count", "GeneRatio", "BgRatio"),
          names(selected)
        )
      ))
      
    } else {
      
      display_cols <- head(names(selected), 6)
    }
    
    print(
      selected[, display_cols, drop = FALSE],
      row.names = FALSE
    )
    
    selected$source_file <- filename
    
    pathway_results[[filename]] <- selected
    
  } else {
    
    cat("Nessuna colonna statistica riconosciuta.\n")
    print(head(dat, 5))
  }
}

# ------------------------------------------------------------
# 3. ESPORTAZIONE RISULTATI PRINCIPALI
# ------------------------------------------------------------

if (length(pathway_results) > 0) {
  
  combined_results <- do.call(
    rbind,
    lapply(pathway_results, function(x) {
      x
    })
  )
  
  write.csv(
    combined_results,
    file.path(output_dir, "R27_top_enrichment_results.csv"),
    row.names = FALSE
  )
  
  cat("\nRisultati estratti:", nrow(combined_results), "\n")
}

# ------------------------------------------------------------
# 4. RIEPILOGO ENRICHMENT
# ------------------------------------------------------------

cat("\n4. RIEPILOGO ENRICHMENT\n")
cat("------------------------------------------------------------\n")

summary_file <- file.path(
  base_dir,
  "functional_enrichment_summary_strict.csv"
)

if (file.exists(summary_file)) {
  
  enrichment_summary <- read.csv(
    summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("File:", basename(summary_file), "\n")
  cat("Dimensioni:", nrow(enrichment_summary),
      "righe x", ncol(enrichment_summary), "colonne\n")
  
  print(
    enrichment_summary,
    row.names = FALSE
  )
  
} else {
  
  cat("Riepilogo strict non trovato.\n")
}

# ------------------------------------------------------------
# 5. INVENTARIO FIGURE
# ------------------------------------------------------------

cat("\n5. FIGURE ESISTENTI\n")
cat("------------------------------------------------------------\n")

if (dir.exists(plots_dir)) {
  
  plot_files <- list.files(
    plots_dir,
    pattern = "\\.(png|pdf|svg|jpg|jpeg|tiff)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(plot_files) > 0) {
    
    plot_inventory <- data.frame(
      file = basename(plot_files),
      relative_path = substring(
        plot_files,
        nchar(plots_dir) + 2
      ),
      extension = tools::file_ext(plot_files),
      size_kb = round(
        file.info(plot_files)$size / 1024,
        2
      ),
      stringsAsFactors = FALSE
    )
    
    print(plot_inventory, row.names = FALSE)
    
    write.csv(
      plot_inventory,
      file.path(output_dir, "R27_existing_figures.csv"),
      row.names = FALSE
    )
    
  } else {
    
    cat("Nessuna figura trovata.\n")
  }
  
} else {
  
  cat("Cartella plots non trovata.\n")
}

# ------------------------------------------------------------
# 6. REPORT DI ROBUSTEZZA
# ------------------------------------------------------------

cat("\n6. REPORT DI ROBUSTEZZA\n")
cat("------------------------------------------------------------\n")

robustness_files <- list.files(
  reports_dir,
  pattern = "\\.(csv|md)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

robustness_files <- robustness_files[
  grepl(
    "robustness|pathway|mapping",
    basename(robustness_files),
    ignore.case = TRUE
  )
]

if (length(robustness_files) > 0) {
  
  cat(
    paste(basename(robustness_files), collapse = "\n"),
    "\n"
  )
  
} else {
  
  cat("Nessun report di robustezza aggiuntivo trovato.\n")
}

# ------------------------------------------------------------
# 7. COMPLETAMENTO
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("R27 INVENTORY COMPLETED\n")
cat("============================================================\n")

cat("\nFile di enrichment:", length(result_files), "\n")
cat("Figure individuate:",
    if (exists("plot_files")) length(plot_files) else 0,
    "\n")

cat("\nOutput inventario:\n")
cat(normalizePath(output_dir, winslash = "/"), "\n")

cat("\nIncolla qui l'intero output della console.\n")
