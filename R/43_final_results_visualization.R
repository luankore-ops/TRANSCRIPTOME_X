
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 43 - FINAL RESULTS VISUALIZATION
# CORRECTED VERSION
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X - SCRIPT 43\n")
cat("FINAL RESULTS VISUALIZATION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(project_dir, "results")

input_dir <- file.path(
  results_dir,
  "final_biological_interpretation",
  "tables"
)

network_dir <- file.path(
  results_dir,
  "gene_pathway_module_network",
  "tables"
)

output_dir <- file.path(
  results_dir,
  "final_results_visualization"
)

png_dir <- file.path(output_dir, "png")
pdf_dir <- file.path(output_dir, "pdf")
audit_dir <- file.path(output_dir, "audit")

for (d in c(png_dir, pdf_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ranking_file <- file.path(
  input_dir, "final_integrated_candidate_ranking.csv"
)

shortlist_file <- file.path(
  input_dir, "final_candidate_shortlist.csv"
)

robust_file <- file.path(
  input_dir, "robust_cross_layer_candidate_profile.csv"
)

evidence_file <- file.path(
  input_dir, "final_candidate_evidence_matrix.csv"
)


# ------------------------------------------------------------
# 2. FUNCTIONS
# ------------------------------------------------------------

read_csv <- function(path) {
  
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN", "N/A")
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

check_columns <- function(data, required, label) {
  
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(
      "Missing columns in ", label, ": ",
      paste(missing, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(data), collapse = ", ")
    )
  }
}

find_column <- function(data, candidates) {
  
  normalized <- tolower(
    gsub("[^a-z0-9]", "", names(data))
  )
  
  targets <- tolower(
    gsub("[^a-z0-9]", "", candidates)
  )
  
  idx <- match(targets, normalized)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) return(NULL)
  
  names(data)[idx[1]]
}

clean_boolean <- function(x) {
  
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  
  tolower(trimws(as.character(x))) %in%
    c("true", "t", "1", "yes")
}

save_figure <- function(name, plot_function,
                        width = 10, height = 7) {
  
  png_path <- file.path(png_dir, paste0(name, ".png"))
  pdf_path <- file.path(pdf_dir, paste0(name, ".pdf"))
  
  grDevices::png(
    filename = png_path,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  
  tryCatch(
    plot_function(),
    finally = grDevices::dev.off()
  )
  
  grDevices::pdf(
    file = pdf_path,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  tryCatch(
    plot_function(),
    finally = grDevices::dev.off()
  )
  
  if (!file.exists(png_path) || !file.exists(pdf_path)) {
    stop("Figure export failed: ", name)
  }
  
  cat("Generated:", name, "\n")
}


# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------

cat("Loading Script 42 results...\n")

ranking <- read_csv(ranking_file)
shortlist <- read_csv(shortlist_file)
robust <- read_csv(robust_file)
evidence <- read_csv(evidence_file)

cat("Integrated candidates:", nrow(ranking), "\n")
cat("Shortlist:", nrow(shortlist), "\n")
cat("Robust candidates:", nrow(robust), "\n\n")


# ------------------------------------------------------------
# 4. INPUT VALIDATION
# ------------------------------------------------------------

ranking_required <- c(
  "gene_symbol",
  "validation_score",
  "validation_class",
  "logFC_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_3",
  "FDR_siNRP1_5"
)

shortlist_required <- ranking_required

# IMPORTANT:
# Script 42 robust profile does not contain validation_class.
# All rows in this file are already ROBUST_CROSS_LAYER.

robust_required <- c(
  "gene_symbol",
  "validation_score",
  "logFC_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_3",
  "FDR_siNRP1_5",
  "expression_pattern",
  "pathway_count",
  "module_count"
)

evidence_required <- c(
  "gene_symbol",
  "network_evidence",
  "pathway_evidence",
  "module_evidence"
)

check_columns(ranking, ranking_required, "ranking")
check_columns(shortlist, shortlist_required, "shortlist")
check_columns(robust, robust_required, "robust profile")
check_columns(evidence, evidence_required, "evidence matrix")

for (col in c(
  "validation_score",
  "logFC_siNRP1_3",
  "logFC_siNRP1_5",
  "FDR_siNRP1_3",
  "FDR_siNRP1_5"
)) {
  ranking[[col]] <- safe_numeric(ranking[[col]])
  shortlist[[col]] <- safe_numeric(shortlist[[col]])
  robust[[col]] <- safe_numeric(robust[[col]])
}

if (anyDuplicated(ranking$gene_symbol) ||
    anyDuplicated(shortlist$gene_symbol) ||
    anyDuplicated(robust$gene_symbol)) {
  stop("Duplicated gene symbols detected.")
}

if (any(!is.finite(robust$validation_score))) {
  stop("Invalid robust candidate scores.")
}

if (any(!is.finite(robust$logFC_siNRP1_3)) ||
    any(!is.finite(robust$logFC_siNRP1_5))) {
  stop("Invalid robust candidate logFC values.")
}

if (any(
  !is.na(robust$FDR_siNRP1_3) &
  (robust$FDR_siNRP1_3 < 0 |
   robust$FDR_siNRP1_3 > 1)
) || any(
  !is.na(robust$FDR_siNRP1_5) &
  (robust$FDR_siNRP1_5 < 0 |
   robust$FDR_siNRP1_5 > 1)
)) {
  stop("Invalid FDR values in robust profile.")
}

if (nrow(robust) != 23 || nrow(shortlist) != 79) {
  warning("Candidate counts differ from previous Script 42 results.")
}

cat("Input validation passed.\n\n")


# ------------------------------------------------------------
# 5. FIGURE 1 - ROBUST CANDIDATE RANKING
# ------------------------------------------------------------

cat("Generating Figure 1...\n")

robust <- robust[
  order(robust$validation_score),
  ,
  drop = FALSE
]

save_figure(
  "Figure_01_Robust_Candidate_Ranking",
  function() {
    
    par(mar = c(5, 9, 4, 2), las = 1)
    
    scores <- robust$validation_score
    labels <- robust$gene_symbol
    
    positions <- barplot(
      scores,
      names.arg = labels,
      horiz = TRUE,
      col = "#3977A8",
      border = NA,
      xlim = c(0, max(scores) + 12),
      xlab = "Integrated validation score",
      main = "Robust cross-layer candidate ranking",
      cex.names = 0.85
    )
    
    text(
      scores + 1,
      positions,
      labels = sprintf("%.1f", scores),
      pos = 4,
      cex = 0.72
    )
  },
  width = 10,
  height = 9
)


# ------------------------------------------------------------
# 6. FIGURE 2 - EXPRESSION HEATMAP
# ------------------------------------------------------------

cat("Generating Figure 2...\n")

heatmap_data <- shortlist[
  order(-shortlist$validation_score),
  ,
  drop = FALSE
]

expression_matrix <- as.matrix(
  heatmap_data[
    ,
    c("logFC_siNRP1_3", "logFC_siNRP1_5")
  ]
)

rownames(expression_matrix) <- heatmap_data$gene_symbol

heatmap_limit <- max(abs(expression_matrix), na.rm = TRUE)

if (!is.finite(heatmap_limit) || heatmap_limit == 0) {
  heatmap_limit <- 1
}

heatmap_limit <- min(heatmap_limit, 5)

save_figure(
  "Figure_02_Shortlist_Expression_Heatmap",
  function() {
    
    n <- nrow(expression_matrix)
    
    par(mar = c(7, 9, 4, 3), xpd = NA)
    
    plot(
      NA,
      xlim = c(0.5, 2.5),
      ylim = c(0.5, n + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "Differential expression across NRP1 knockdowns"
    )
    
    axis(
      1,
      at = 1:2,
      labels = c("siNRP1.3", "siNRP1.5")
    )
    
    axis(
      2,
      at = seq_len(n),
      labels = rev(rownames(expression_matrix)),
      las = 1,
      cex.axis = 0.48
    )
    
    palette <- colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(101)
    
    for (i in seq_len(n)) {
      
      row_index <- n - i + 1
      
      for (j in 1:2) {
        
        value <- expression_matrix[i, j]
        
        if (is.na(value)) {
          fill <- "grey85"
        } else {
          
          scaled <- max(-1, min(1, value / heatmap_limit))
          
          color_index <- round((scaled + 1) * 50) + 1
          
          fill <- palette[color_index]
        }
        
        rect(
          j - 0.42,
          row_index - 0.42,
          j + 0.42,
          row_index + 0.42,
          col = fill,
          border = "white"
        )
        
        fdr <- if (j == 1) {
          heatmap_data$FDR_siNRP1_3[i]
        } else {
          heatmap_data$FDR_siNRP1_5[i]
        }
        
        if (!is.na(fdr) && fdr < 0.05) {
          points(
            j,
            row_index,
            pch = 16,
            cex = 0.55,
            col = "black"
          )
        }
      }
    }
    
    mtext(
      "Color = logFC; black dot = FDR < 0.05",
      side = 1,
      line = 4.5,
      cex = 0.8
    )
    
    legend(
      "topright",
      legend = c("Downregulated", "No change", "Upregulated"),
      fill = c("#2166AC", "white", "#B2182B"),
      border = "grey60",
      cex = 0.7,
      bg = "white"
    )
  },
  width = 9,
  height = 14
)


# ------------------------------------------------------------
# 7. FIGURE 3 - INTEGRATED EVIDENCE MATRIX
# ------------------------------------------------------------

cat("Generating Figure 3...\n")

top_genes <- head(
  shortlist$gene_symbol[
    order(-shortlist$validation_score)
  ],
  25
)

top_evidence <- evidence[
  match(top_genes, evidence$gene_symbol),
  ,
  drop = FALSE
]

top_evidence <- top_evidence[
  !is.na(top_evidence$gene_symbol),
  ,
  drop = FALSE
]

evidence_matrix <- sapply(
  c(
    "network_evidence",
    "pathway_evidence",
    "module_evidence"
  ),
  function(x) clean_boolean(top_evidence[[x]])
)

rownames(evidence_matrix) <- top_evidence$gene_symbol

save_figure(
  "Figure_03_Integrated_Evidence_Matrix",
  function() {
    
    n <- nrow(evidence_matrix)
    
    par(mar = c(6, 8, 4, 2), xpd = NA)
    
    plot(
      NA,
      xlim = c(0.5, 3.5),
      ylim = c(0.5, n + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "Integrated evidence for top candidates"
    )
    
    axis(
      1,
      at = 1:3,
      labels = c("Network", "Pathway", "Module")
    )
    
    axis(
      2,
      at = seq_len(n),
      labels = rev(rownames(evidence_matrix)),
      las = 1,
      cex.axis = 0.8
    )
    
    for (i in seq_len(n)) {
      
      row_index <- n - i + 1
      
      for (j in 1:3) {
        
        present <- evidence_matrix[i, j]
        
        rect(
          j - 0.38,
          row_index - 0.38,
          j + 0.38,
          row_index + 0.38,
          col = if (present) "#3977A8" else "grey90",
          border = "white"
        )
        
        text(
          j,
          row_index,
          labels = if (present) "+" else "-",
          col = if (present) "white" else "grey40",
          cex = 0.9
        )
      }
    }
  },
  width = 9,
  height = 9
)


# ------------------------------------------------------------
# 8. FIGURE 4 - GENE-MODULE NETWORK
# ------------------------------------------------------------

cat("Generating Figure 4...\n")

gene_module_file <- file.path(
  network_dir,
  "gene_module_edges.csv"
)

network_generated <- FALSE

if (file.exists(gene_module_file)) {
  
  edges <- read_csv(gene_module_file)
  
  gene_col <- find_column(
    edges,
    c("gene_symbol", "gene", "symbol", "gene_name")
  )
  
  module_col <- find_column(
    edges,
    c("module", "module_name", "module_id", "gene_module")
  )
  
  if (!is.null(gene_col) && !is.null(module_col)) {
    
    edges[[gene_col]] <- as.character(edges[[gene_col]])
    edges[[module_col]] <- as.character(edges[[module_col]])
    
    edges <- edges[
      !is.na(edges[[gene_col]]) &
        !is.na(edges[[module_col]]) &
        edges[[gene_col]] != "" &
        edges[[module_col]] != "",
      ,
      drop = FALSE
    ]
    
    edges <- unique(
      edges[
        edges[[gene_col]] %in% robust$gene_symbol,
        ,
        drop = FALSE
      ]
    )
    
    if (nrow(edges) > 0) {
      
      genes <- sort(unique(edges[[gene_col]]))
      modules <- sort(unique(edges[[module_col]]))
      
      gene_y <- seq(
        length(genes),
        1,
        length.out = length(genes)
      )
      
      module_y <- seq(
        length(modules),
        1,
        length.out = length(modules)
      )
      
      gene_positions <- setNames(gene_y, genes)
      module_positions <- setNames(module_y, modules)
      
      save_figure(
        "Figure_04_Robust_Gene_Module_Network",
        function() {
          
          par(mar = c(3, 2, 4, 2), xpd = NA)
          
          plot(
            NA,
            xlim = c(0, 10),
            ylim = c(
              0,
              max(length(genes), length(modules)) + 1
            ),
            xaxt = "n",
            yaxt = "n",
            xlab = "",
            ylab = "",
            main = "Robust candidates and gene-module relationships"
          )
          
          for (i in seq_len(nrow(edges))) {
            
            gene <- edges[[gene_col]][i]
            module <- edges[[module_col]][i]
            
            segments(
              3,
              gene_positions[gene],
              7,
              module_positions[module],
              col = "grey75",
              lwd = 1
            )
          }
          
          points(
            rep(3, length(genes)),
            gene_y,
            pch = 21,
            bg = "#3977A8",
            col = "white",
            cex = 1.5
          )
          
          text(
            rep(2.7, length(genes)),
            gene_y,
            labels = genes,
            pos = 2,
            cex = 0.75
          )
          
          points(
            rep(7, length(modules)),
            module_y,
            pch = 22,
            bg = "#D98C3F",
            col = "white",
            cex = 2
          )
          
          text(
            rep(7.3, length(modules)),
            module_y,
            labels = modules,
            pos = 4,
            cex = 0.75
          )
          
          text(
            3,
            max(gene_y) + 0.8,
            "ROBUST GENES",
            font = 2,
            cex = 0.85
          )
          
          text(
            7,
            max(module_y) + 0.8,
            "BIOLOGICAL MODULES",
            font = 2,
            cex = 0.85
          )
        },
        width = 12,
        height = max(7, length(genes) * 0.32)
      )
      
      network_generated <- TRUE
    }
  }
}

if (!network_generated) {
  cat(
    "Figure 4 skipped: compatible gene-module edges not found.\n"
  )
}


# ------------------------------------------------------------
# 9. FIGURE 5 - CANDIDATE CLASS DISTRIBUTION
# ------------------------------------------------------------

cat("Generating Figure 5...\n")

class_counts <- table(
  factor(
    ranking$validation_class,
    levels = c(
      "ROBUST_CROSS_LAYER",
      "MULTI_LAYER_SUPPORTED",
      "TRANSCRIPTOMICALLY_SUPPORTED",
      "NETWORK_PRIORITY_ONLY",
      "EXPLORATORY"
    )
  )
)

class_counts <- class_counts[class_counts > 0]

save_figure(
  "Figure_05_Integrated_Candidate_Class_Distribution",
  function() {
    
    par(mar = c(5, 12, 4, 2), las = 1)
    
    values <- as.numeric(class_counts)
    
    positions <- barplot(
      values,
      names.arg = names(class_counts),
      horiz = TRUE,
      col = c(
        "#245A83",
        "#3977A8",
        "#79A6C7",
        "#B9CBD9",
        "#D9D9D9"
      )[seq_along(values)],
      border = NA,
      xlab = "Number of genes",
      main = "Integrated candidate classification",
      cex.names = 0.75
    )
    
    text(
      values + max(values) * 0.015,
      positions,
      labels = values,
      pos = 4,
      cex = 0.8
    )
  },
  width = 11,
  height = 7
)


# ------------------------------------------------------------
# 10. FIGURE 6 - KNOCKDOWN EFFECT COMPARISON
# ------------------------------------------------------------

cat("Generating Figure 6...\n")

save_figure(
  "Figure_06_Knockdown_LogFC_Comparison",
  function() {
    
    x <- robust$logFC_siNRP1_3
    y <- robust$logFC_siNRP1_5
    
    valid <- is.finite(x) & is.finite(y)
    
    x_valid <- x[valid]
    y_valid <- y[valid]
    labels <- robust$gene_symbol[valid]
    
    plot(
      x_valid,
      y_valid,
      pch = 21,
      bg = "#3977A8",
      col = "white",
      cex = 1.5,
      xlab = "logFC siNRP1.3 vs siNT",
      ylab = "logFC siNRP1.5 vs siNT",
      main = "Expression effect comparison across knockdowns"
    )
    
    abline(
      h = 0,
      v = 0,
      lty = 2,
      col = "grey55"
    )
    
    abline(
      a = 0,
      b = 1,
      lty = 3,
      col = "grey65"
    )
    
    text(
      x_valid,
      y_valid,
      labels = labels,
      pos = 3,
      cex = 0.65
    )
    
    correlation <- if (length(x_valid) >= 3) {
      cor(x_valid, y_valid, use = "complete.obs")
    } else {
      NA_real_
    }
    
    legend(
      "topleft",
      legend = paste(
        "Pearson r =",
        ifelse(is.finite(correlation),
               round(correlation, 3), "NA")
      ),
      bty = "n",
      cex = 0.8
    )
  },
  width = 9,
  height = 8
)


# ------------------------------------------------------------
# 11. FIGURE 7 - TOP 15 CANDIDATES
# ------------------------------------------------------------

cat("Generating Figure 7...\n")

top_candidates <- shortlist[
  order(-shortlist$validation_score),
  ,
  drop = FALSE
]

top_candidates <- head(top_candidates, 15)

top_candidates <- top_candidates[
  order(top_candidates$validation_score),
  ,
  drop = FALSE
]

save_figure(
  "Figure_07_Top_15_Integrated_Candidates",
  function() {
    
    par(mar = c(5, 8, 4, 2), las = 1)
    
    positions <- barplot(
      top_candidates$validation_score,
      names.arg = top_candidates$gene_symbol,
      horiz = TRUE,
      col = "#3977A8",
      border = NA,
      xlab = "Integrated validation score",
      main = "Top 15 integrated candidates",
      xlim = c(
        0,
        max(top_candidates$validation_score) + 12
      ),
      cex.names = 0.8
    )
    
    text(
      top_candidates$validation_score + 1,
      positions,
      labels = round(
        top_candidates$validation_score,
        1
      ),
      pos = 4,
      cex = 0.75
    )
  },
  width = 10,
  height = 8
)


# ------------------------------------------------------------
# 12. FIGURE INVENTORY
# ------------------------------------------------------------

png_files <- list.files(
  png_dir,
  pattern = "\\.png$",
  full.names = FALSE
)

pdf_files <- list.files(
  pdf_dir,
  pattern = "\\.pdf$",
  full.names = FALSE
)

figure_names <- sub("\\.png$", "", png_files)

inventory <- data.frame(
  figure = figure_names,
  png_available = TRUE,
  pdf_available = paste0(figure_names, ".pdf") %in% pdf_files,
  stringsAsFactors = FALSE
)

write.csv(
  inventory,
  file.path(audit_dir, "script43_figure_inventory.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. FINAL AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "integrated_candidates",
    "final_shortlist",
    "robust_candidates",
    "png_figures",
    "pdf_figures",
    "gene_module_network_generated"
  ),
  value = c(
    nrow(ranking),
    nrow(shortlist),
    nrow(robust),
    length(png_files),
    length(pdf_files),
    network_generated
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(audit_dir, "script43_audit.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. COMPLETION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 43 COMPLETED\n")
cat("============================================================\n\n")

cat("Integrated candidates:", nrow(ranking), "\n")
cat("Final shortlist:", nrow(shortlist), "\n")
cat("Robust candidates:", nrow(robust), "\n")
cat("PNG figures:", length(png_files), "\n")
cat("PDF figures:", length(pdf_files), "\n")
cat("Gene-module network generated:", network_generated, "\n")

cat("\nPNG DIRECTORY:\n", png_dir, "\n")
cat("\nPDF DIRECTORY:\n", pdf_dir, "\n")

cat("\n============================================================\n")
