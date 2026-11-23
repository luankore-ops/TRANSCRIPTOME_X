
# ============================================================
# TRANSCRIPTOME-X
# STEP 5 - VOLCANO PLOTS
# ============================================================

cat("\nTRANSCRIPTOME-X - VOLCANO PLOTS\n")
cat("========================================\n")

# 1. PACCHETTI
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

# 2. PERCORSI
# ------------------------------------------------------------

input_dir <- file.path(
  "data",
  "processed",
  "differential_expression",
  "annotated"
)

output_dir <- file.path(
  "data",
  "processed",
  "visualizations",
  "volcano_plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. FUNZIONE PER GENERARE VOLCANO PLOT
# ------------------------------------------------------------

create_volcano <- function(comparison, input_file) {
  
  cat("\nConfronto:", comparison, "\n")
  
  if (!file.exists(input_file)) {
    stop("File non trovato: ", input_file)
  }
  
  data <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )
  
  required_columns <- c(
    "gene_id",
    "gene_name",
    "logFC",
    "FDR"
  )
  
  if (!all(required_columns %in% colnames(data))) {
    stop("Colonne richieste mancanti.")
  }
  
  if (any(!is.finite(data$logFC)) ||
      any(!is.finite(data$FDR))) {
    stop("Sono presenti valori non finiti.")
  }
  
  if (any(data$FDR < 0 | data$FDR > 1)) {
    stop("Valori FDR non validi.")
  }
  
  # Classificazione dei geni
  data$category <- "Not significant"
  
  significant <- data$FDR < 0.05
  
  data$category[
    significant & abs(data$logFC) < 1
  ] <- "Significant"
  
  data$category[
    significant & data$logFC >= 1
  ] <- "Up-regulated"
  
  data$category[
    significant & data$logFC <= -1
  ] <- "Down-regulated"
  
  data$category <- factor(
    data$category,
    levels = c(
      "Not significant",
      "Significant",
      "Up-regulated",
      "Down-regulated"
    )
  )
  
  # Calcolo asse Y
  data$minus_log10_FDR <- -log10(
    pmax(data$FDR, .Machine$double.xmin)
  )
  
  # Conteggi
  category_counts <- as.data.frame(
    table(data$category)
  )
  
  cat("\nDistribuzione categorie:\n")
  print(category_counts)
  
  # Evidenziazione NRP1
  nrp1 <- data[
    data$gene_name == "NRP1",
  ]
  
  if (nrow(nrp1) != 1) {
    stop("NRP1 non identificato univocamente.")
  }
  
  cat("\nNRP1:\n")
  print(
    nrp1[
      ,
      c("gene_name", "logFC", "FDR")
    ]
  )
  
  # Grafico
  volcano <- ggplot(
    data,
    aes(
      x = logFC,
      y = minus_log10_FDR,
      color = category
    )
  ) +
    
    geom_point(
      alpha = 0.65,
      size = 1.2
    ) +
    
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      color = "grey40"
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey40"
    ) +
    
    geom_point(
      data = nrp1,
      aes(x = logFC, y = minus_log10_FDR),
      inherit.aes = FALSE,
      color = "black",
      fill = "gold",
      shape = 23,
      size = 4,
      stroke = 1
    ) +
    
    geom_text(
      data = nrp1,
      aes(
        x = logFC,
        y = minus_log10_FDR,
        label = gene_name
      ),
      inherit.aes = FALSE,
      nudge_y = 0.5,
      fontface = "bold",
      size = 4
    ) +
    
    scale_color_manual(
      values = c(
        "Not significant" = "grey75",
        "Significant" = "steelblue",
        "Up-regulated" = "firebrick2",
        "Down-regulated" = "royalblue3"
      )
    ) +
    
    labs(
      title = paste(
        "Differential Expression:",
        comparison
      ),
      subtitle = paste(
        "edgeR Quasi-Likelihood |",
        "FDR < 0.05 | |log2FC| >= 1"
      ),
      x = "log2 Fold Change",
      y = expression(-log[10]("FDR")),
      color = "Gene category"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10
      ),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # Salvataggio PNG
  ggsave(
    filename = file.path(
      output_dir,
      paste0(comparison, "_volcano.png")
    ),
    plot = volcano,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  # Salvataggio PDF
  ggsave(
    filename = file.path(
      output_dir,
      paste0(comparison, "_volcano.pdf")
    ),
    plot = volcano,
    width = 10,
    height = 7
  )
  
  # Salvataggio dati classificati
  write.csv(
    data,
    file.path(
      output_dir,
      paste0(comparison, "_volcano_data.csv")
    ),
    row.names = FALSE
  )
  
  cat("Grafici PNG e PDF salvati.\n")
  
  return(
    data.frame(
      comparison = comparison,
      genes_tested = nrow(data),
      not_significant = sum(
        data$category == "Not significant"
      ),
      significant_below_fc_threshold = sum(
        data$category == "Significant"
      ),
      upregulated = sum(
        data$category == "Up-regulated"
      ),
      downregulated = sum(
        data$category == "Down-regulated"
      ),
      stringsAsFactors = FALSE
    )
  )
}

# 4. GENERAZIONE DEI DUE GRAFICI
# ------------------------------------------------------------

summary_3 <- create_volcano(
  "siNRP1_3_vs_siNT",
  file.path(
    input_dir,
    "siNRP1_3_vs_siNT_annotated.csv"
  )
)

summary_5 <- create_volcano(
  "siNRP1_5_vs_siNT",
  file.path(
    input_dir,
    "siNRP1_5_vs_siNT_annotated.csv"
  )
)

# 5. RIEPILOGO
# ------------------------------------------------------------

summary_results <- rbind(
  summary_3,
  summary_5
)

write.csv(
  summary_results,
  file.path(
    output_dir,
    "volcano_summary.csv"
  ),
  row.names = FALSE
)

cat("\nRIEPILOGO VOLCANO PLOTS\n")
cat("----------------------------\n")

print(summary_results)

# 6. VALIDAZIONE OUTPUT
# ------------------------------------------------------------

expected_files <- c(
  "siNRP1_3_vs_siNT_volcano.png",
  "siNRP1_3_vs_siNT_volcano.pdf",
  "siNRP1_3_vs_siNT_volcano_data.csv",
  "siNRP1_5_vs_siNT_volcano.png",
  "siNRP1_5_vs_siNT_volcano.pdf",
  "siNRP1_5_vs_siNT_volcano_data.csv",
  "volcano_summary.csv"
)

file_status <- file.exists(
  file.path(output_dir, expected_files)
)

cat("\nCONTROLLO FILE\n")
print(
  data.frame(
    file = expected_files,
    exists = file_status
  )
)

if (!all(file_status)) {
  stop("Uno o più file non sono stati generati.")
}

cat("\n========================================\n")
cat("STEP 5 COMPLETATO\n")
cat("Volcano plots generati correttamente.\n")
cat("========================================\n")
