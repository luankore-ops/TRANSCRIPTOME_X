
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 39 - NETWORK VISUALIZATION AND BIOLOGICAL INTERPRETATION
#
# Compatible with Script 37 and Script 38 output schemas.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

cat("\n============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 39 - NETWORK VISUALIZATION AND BIOLOGICAL INTERPRETATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------

project_dir <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = TRUE
)

network_dir <- file.path(
  project_dir, "results", "gene_pathway_module_network"
)

topology_dir <- file.path(
  project_dir, "results", "network_topology_analysis"
)

output_dir <- file.path(
  project_dir, "results", "network_visualization_analysis"
)

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
report_dir <- file.path(output_dir, "reports")
audit_dir <- file.path(output_dir, "audit")

for (d in c(table_dir, figure_dir, report_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# 2. DEPENDENCIES
# ------------------------------------------------------------

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required.")
}

library(igraph)

# ------------------------------------------------------------
# 3. INPUT FILES
# ------------------------------------------------------------

files <- c(
  gene_pathway = file.path(
    network_dir, "tables", "gene_pathway_edges.csv"
  ),
  gene_module = file.path(
    network_dir, "tables", "gene_module_edges.csv"
  ),
  pathway_module = file.path(
    network_dir, "tables", "pathway_module_edges.csv"
  ),
  node_topology = file.path(
    topology_dir, "tables", "node_topology_metrics.csv"
  ),
  gene_ranking = file.path(
    topology_dir, "tables", "gene_topology_ranking.csv"
  ),
  bridge_topology = file.path(
    topology_dir, "tables", "bridge_gene_topology.csv"
  ),
  community_assignments = file.path(
    topology_dir, "tables", "node_community_assignments.csv"
  ),
  community_summary = file.path(
    topology_dir, "tables", "network_community_summary.csv"
  ),
  component_summary = file.path(
    topology_dir, "tables", "connected_component_summary.csv"
  )
)

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

read_csv <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

gene_pathway <- read_csv(files["gene_pathway"])
gene_module <- read_csv(files["gene_module"])
pathway_module <- read_csv(files["pathway_module"])

node_topology <- read_csv(files["node_topology"])
gene_ranking <- read_csv(files["gene_ranking"])
bridge_topology <- read_csv(files["bridge_topology"])
community_assignments <- read_csv(files["community_assignments"])
community_summary <- read_csv(files["community_summary"])
component_summary <- read_csv(files["component_summary"])

# ------------------------------------------------------------
# 4. VALIDATE ACTUAL SCRIPT 37 SCHEMAS
# ------------------------------------------------------------

required_columns <- list(
  gene_pathway = c("pathway_id", "gene_symbol"),
  gene_module = c("gene_symbol", "module"),
  pathway_module = c("pathway_id", "module")
)

for (nm in names(required_columns)) {
  
  dat <- get(nm)
  missing <- setdiff(required_columns[[nm]], colnames(dat))
  
  if (length(missing) > 0) {
    stop(
      nm, " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

required_topology <- c(
  "node_id", "node_type", "label", "degree",
  "degree_centrality", "betweenness",
  "harmonic_centrality", "component_id",
  "community_id", "is_articulation_point"
)

missing_topology <- setdiff(
  required_topology,
  colnames(node_topology)
)

if (length(missing_topology) > 0) {
  stop(
    "node_topology_metrics.csv is missing: ",
    paste(missing_topology, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. BUILD EDGES USING REAL COLUMN NAMES
# ------------------------------------------------------------

gene_pathway_edges <- data.frame(
  from = paste0("GENE::", gene_pathway$gene_symbol),
  to = paste0("PATHWAY::", gene_pathway$pathway_id),
  edge_type = "GENE_PATHWAY"
)

gene_module_edges <- data.frame(
  from = paste0("GENE::", gene_module$gene_symbol),
  to = paste0("MODULE::", gene_module$module),
  edge_type = "GENE_MODULE"
)

pathway_module_edges <- data.frame(
  from = paste0("PATHWAY::", pathway_module$pathway_id),
  to = paste0("MODULE::", pathway_module$module),
  edge_type = "PATHWAY_MODULE"
)

edges <- rbind(
  gene_pathway_edges,
  gene_module_edges,
  pathway_module_edges
)

edges$from <- trimws(as.character(edges$from))
edges$to <- trimws(as.character(edges$to))

edges <- edges[
  !is.na(edges$from) &
    !is.na(edges$to) &
    nzchar(edges$from) &
    nzchar(edges$to) &
    edges$from != edges$to,
  ,
  drop = FALSE
]

edges <- unique(edges)

if (nrow(edges) == 0) {
  stop("No valid network edges were generated.")
}

# ------------------------------------------------------------
# 6. BUILD NODE TABLE
# ------------------------------------------------------------

nodes <- unique(c(edges$from, edges$to))

node_table <- data.frame(
  name = nodes,
  node_type = sub("::.*$", "", nodes),
  label = sub("^[^:]+::", "", nodes),
  stringsAsFactors = FALSE
)

# Map topology using node type and label.
topology_map <- node_topology

topology_map$node_type <- toupper(
  as.character(topology_map$node_type)
)

topology_map$label <- as.character(topology_map$label)

node_table$topology_id <- paste0(
  node_table$node_type,
  "::",
  node_table$label
)

topology_map$topology_id <- paste0(
  topology_map$node_type,
  "::",
  topology_map$label
)

topology_columns <- c(
  "topology_id",
  "node_id",
  "degree",
  "degree_centrality",
  "betweenness",
  "harmonic_centrality",
  "component_id",
  "community_id",
  "is_articulation_point"
)

topology_columns <- intersect(
  topology_columns,
  colnames(topology_map)
)

node_table <- merge(
  node_table,
  topology_map[, topology_columns, drop = FALSE],
  by = "topology_id",
  all.x = TRUE,
  sort = FALSE
)

node_table <- node_table[
  match(nodes, node_table$name),
  ,
  drop = FALSE
]

if (anyDuplicated(node_table$name)) {
  stop("Duplicate node identifiers detected.")
}

# ------------------------------------------------------------
# 7. BUILD GRAPH
# ------------------------------------------------------------

graph <- graph_from_data_frame(
  d = edges[, c("from", "to", "edge_type")],
  directed = FALSE,
  vertices = node_table
)

cat("Network reconstructed successfully.\n")
cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n\n")

# ------------------------------------------------------------
# 8. EXPORT NODE AND EDGE TABLES
# ------------------------------------------------------------

node_export <- as_data_frame(graph, what = "vertices")
edge_export <- as_data_frame(graph, what = "edges")

write.csv(
  node_export,
  file.path(table_dir, "network_visualization_nodes.csv"),
  row.names = FALSE
)

write.csv(
  edge_export,
  file.path(table_dir, "network_visualization_edges.csv"),
  row.names = FALSE
)

edge_type_summary <- as.data.frame(
  table(edges$edge_type),
  stringsAsFactors = FALSE
)

colnames(edge_type_summary) <- c(
  "edge_type", "edge_count"
)

edge_type_summary <- edge_type_summary[
  order(-edge_type_summary$edge_count),
]

write.csv(
  edge_type_summary,
  file.path(table_dir, "network_edge_type_composition.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. CENTRALITY RANKING
# ------------------------------------------------------------

centrality_table <- node_export

centrality_table <- centrality_table[
  order(
    -centrality_table$betweenness,
    -centrality_table$degree,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

centrality_table$betweenness_rank <- seq_len(
  nrow(centrality_table)
)

write.csv(
  centrality_table,
  file.path(table_dir, "network_centrality_ranking.csv"),
  row.names = FALSE
)

top_30 <- head(centrality_table, 30)

write.csv(
  top_30,
  file.path(table_dir, "top_30_network_nodes.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. BRIDGE GENES
# ------------------------------------------------------------

if (!"gene_symbol" %in% colnames(bridge_topology)) {
  
  if ("gene" %in% colnames(bridge_topology)) {
    names(bridge_topology)[
      names(bridge_topology) == "gene"
    ] <- "gene_symbol"
  } else {
    stop("Bridge table has no gene_symbol/gene column.")
  }
}

bridge_topology$node_name <- paste0(
  "GENE::",
  bridge_topology$gene_symbol
)

bridge_topology <- bridge_topology[
  bridge_topology$node_name %in% V(graph)$name,
  ,
  drop = FALSE
]

write.csv(
  bridge_topology,
  file.path(table_dir, "bridge_gene_visualization_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. COMPONENT AND COMMUNITY AUDIT
# ------------------------------------------------------------

component_result <- components(graph)

component_table <- data.frame(
  component_id = seq_along(component_result$csize),
  node_count = as.integer(component_result$csize)
)

component_table <- component_table[
  order(-component_table$node_count),
]

write.csv(
  component_table,
  file.path(table_dir, "visualization_component_summary.csv"),
  row.names = FALSE
)

community_table <- node_export[
  !is.na(node_export$community_id),
  c("name", "node_type", "label", "community_id"),
  drop = FALSE
]

community_summary_39 <- as.data.frame(
  table(community_table$community_id),
  stringsAsFactors = FALSE
)

colnames(community_summary_39) <- c(
  "community_id", "node_count"
)

community_summary_39 <- community_summary_39[
  order(-community_summary_39$node_count),
]

write.csv(
  community_table,
  file.path(table_dir, "visualization_community_assignments.csv"),
  row.names = FALSE
)

write.csv(
  community_summary_39,
  file.path(table_dir, "visualization_community_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. FIGURE 1 - FULL NETWORK
# ------------------------------------------------------------

set.seed(39)

network_layout <- layout_with_fr(
  graph,
  niter = 500,
  grid = "nogrid"
)

node_colors <- ifelse(
  V(graph)$node_type == "GENE",
  "#4472A1",
  ifelse(
    V(graph)$node_type == "PATHWAY",
    "#D58A42",
    "#4D967A"
  )
)

node_sizes <- pmin(
  3 + degree(graph) * 0.12,
  9
)

png(
  file.path(figure_dir, "01_full_network_overview.png"),
  width = 2400,
  height = 2000,
  res = 220
)

par(mar = c(1, 1, 4, 1))

plot(
  graph,
  layout = network_layout,
  vertex.color = node_colors,
  vertex.size = node_sizes,
  vertex.frame.color = NA,
  vertex.label = NA,
  edge.color = grDevices::adjustcolor(
    "grey50", alpha.f = 0.16
  ),
  edge.width = 0.35,
  main = "TRANSCRIPTOME-X Integrated Network"
)

legend(
  "topleft",
  legend = c("Genes", "Pathways", "Modules"),
  col = c("#4472A1", "#D58A42", "#4D967A"),
  pch = 19,
  pt.cex = 1.5,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 13. FIGURE 2 - TOP BETWEENNESS
# ------------------------------------------------------------

top_betweenness <- head(
  centrality_table[
    order(-centrality_table$betweenness),
    ,
    drop = FALSE
  ],
  20
)

png(
  file.path(figure_dir, "02_top_betweenness_nodes.png"),
  width = 2000,
  height = 1400,
  res = 200
)

par(mar = c(5, 13, 4, 2))

barplot(
  rev(top_betweenness$betweenness),
  names.arg = rev(top_betweenness$label),
  horiz = TRUE,
  las = 1,
  col = "#4472A1",
  border = NA,
  main = "Top 20 Nodes by Betweenness Centrality",
  xlab = "Betweenness centrality"
)

dev.off()

# ------------------------------------------------------------
# 14. FIGURE 3 - BRIDGE GENES
# ------------------------------------------------------------

if (nrow(bridge_topology) > 0) {
  
  bridge_plot <- bridge_topology
  
  if (!"betweenness" %in% colnames(bridge_plot)) {
    bridge_plot <- merge(
      bridge_plot,
      node_export[, c("name", "betweenness")],
      by.x = "node_name",
      by.y = "name",
      all.x = TRUE
    )
  }
  
  bridge_plot <- bridge_plot[
    order(bridge_plot$betweenness),
    ,
    drop = FALSE
  ]
  
  png(
    file.path(figure_dir, "03_bridge_gene_betweenness.png"),
    width = 1800,
    height = 1100,
    res = 200
  )
  
  par(mar = c(5, 8, 4, 2))
  
  barplot(
    bridge_plot$betweenness,
    names.arg = bridge_plot$gene_symbol,
    horiz = TRUE,
    las = 1,
    col = "#4D967A",
    border = NA,
    main = "Bridge Gene Betweenness Centrality",
    xlab = "Betweenness centrality"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# 15. FIGURE 4 - EDGE COMPOSITION
# ------------------------------------------------------------

png(
  file.path(figure_dir, "04_network_edge_composition.png"),
  width = 1500,
  height = 1200,
  res = 200
)

par(mar = c(8, 5, 4, 2))

barplot(
  edge_type_summary$edge_count,
  names.arg = edge_type_summary$edge_type,
  las = 2,
  col = "#D58A42",
  border = NA,
  main = "Integrated Network Edge Composition",
  ylab = "Number of edges"
)

dev.off()

# ------------------------------------------------------------
# 16. FIGURE 5 - BRIDGE SUBNETWORK
# ------------------------------------------------------------

bridge_names <- intersect(
  bridge_topology$node_name,
  V(graph)$name
)

if (length(bridge_names) > 0) {
  
  bridge_neighbors <- unique(unlist(
    lapply(
      bridge_names,
      function(x) as_ids(neighbors(graph, x))
    )
  ))
  
  focused_nodes <- unique(c(
    bridge_names,
    bridge_neighbors
  ))
  
  focused_graph <- induced_subgraph(
    graph,
    vids = focused_nodes
  )
  
  set.seed(390)
  
  focused_layout <- layout_with_fr(
    focused_graph,
    niter = 800
  )
  
  focused_colors <- ifelse(
    V(focused_graph)$node_type == "GENE",
    "#4472A1",
    ifelse(
      V(focused_graph)$node_type == "PATHWAY",
      "#D58A42",
      "#4D967A"
    )
  )
  
  focused_sizes <- ifelse(
    V(focused_graph)$name %in% bridge_names,
    12,
    5
  )
  
  png(
    file.path(figure_dir, "05_bridge_gene_subnetwork.png"),
    width = 2400,
    height = 2000,
    res = 220
  )
  
  par(mar = c(1, 1, 4, 1))
  
  plot(
    focused_graph,
    layout = focused_layout,
    vertex.color = focused_colors,
    vertex.size = focused_sizes,
    vertex.frame.color = "grey30",
    vertex.label = ifelse(
      V(focused_graph)$name %in% bridge_names,
      V(focused_graph)$label,
      NA
    ),
    vertex.label.cex = 0.9,
    vertex.label.color = "black",
    edge.color = grDevices::adjustcolor(
      "grey40", alpha.f = 0.4
    ),
    edge.width = 0.7,
    main = "Bridge Gene Neighborhood Network"
  )
  
  legend(
    "topleft",
    legend = c("Gene", "Pathway", "Module", "Bridge gene"),
    col = c(
      "#4472A1", "#D58A42", "#4D967A", "grey30"
    ),
    pch = 19,
    pt.cex = c(1.3, 1.3, 1.3, 1.8),
    bty = "n"
  )
  
  dev.off()
  
  write.csv(
    as_data_frame(focused_graph, what = "vertices"),
    file.path(table_dir, "bridge_subnetwork_nodes.csv"),
    row.names = FALSE
  )
  
  write.csv(
    as_data_frame(focused_graph, what = "edges"),
    file.path(table_dir, "bridge_subnetwork_edges.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 17. FIGURE 6 - COMMUNITY DISTRIBUTION
# ------------------------------------------------------------

if (nrow(community_summary_39) > 0) {
  
  community_plot <- head(
    community_summary_39,
    20
  )
  
  png(
    file.path(figure_dir, "06_community_size_distribution.png"),
    width = 1800,
    height = 1200,
    res = 200
  )
  
  par(mar = c(7, 5, 4, 2))
  
  barplot(
    community_plot$node_count,
    names.arg = community_plot$community_id,
    las = 2,
    col = "#7B83A8",
    border = NA,
    main = "Louvain Community Size Distribution",
    ylab = "Number of nodes"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# 18. FIGURE 7 - CONNECTED COMPONENTS
# ------------------------------------------------------------

png(
  file.path(figure_dir, "07_connected_component_sizes.png"),
  width = 1600,
  height = 1100,
  res = 200
)

par(mar = c(5, 5, 4, 2))

barplot(
  component_table$node_count,
  names.arg = component_table$component_id,
  col = "#8C9A83",
  border = NA,
  main = "Connected Component Size Distribution",
  xlab = "Component",
  ylab = "Number of nodes"
)

dev.off()

# ------------------------------------------------------------
# 19. BIOLOGICAL INTERPRETATION REPORT
# ------------------------------------------------------------

top_text <- paste(
  head(
    paste0(
      centrality_table$label,
      " (",
      signif(centrality_table$betweenness, 4),
      ")"
    ),
    10
  ),
  collapse = "\n"
)

bridge_text <- if (nrow(bridge_topology) > 0) {
  paste(
    paste0(
      bridge_topology$gene_symbol,
      " | ",
      if ("module_count" %in% colnames(bridge_topology)) {
        paste0("modules=", bridge_topology$module_count)
      } else {
        "module count unavailable"
      },
      " | betweenness=",
      if ("betweenness" %in% colnames(bridge_topology)) {
        signif(bridge_topology$betweenness, 4)
      } else {
        "NA"
      }
    ),
    collapse = "\n"
  )
} else {
  "No bridge genes available."
}

report <- c(
  "TRANSCRIPTOME-X",
  "NETWORK VISUALIZATION AND BIOLOGICAL INTERPRETATION",
  "===================================================",
  "",
  paste("Network nodes:", vcount(graph)),
  paste("Network edges:", ecount(graph)),
  paste("Connected components:", length(component_result$csize)),
  paste("Communities represented:", nrow(community_summary_39)),
  paste("Bridge genes:", length(bridge_names)),
  "",
  "TOP NODES BY BETWEENNESS CENTRALITY",
  "-----------------------------------",
  top_text,
  "",
  "BRIDGE GENES",
  "------------",
  bridge_text,
  "",
  "INTERPRETATION",
  "--------------",
  "Degree describes the number of direct network connections.",
  "Betweenness describes the extent to which a node lies on shortest paths.",
  "Bridge genes connect multiple annotated modules or network neighborhoods.",
  "A higher module count does not necessarily imply higher betweenness.",
  "Network communities are computational partitions, not necessarily",
  "independent biological pathways.",
  "",
  "IMPORTANT LIMITATIONS",
  "---------------------",
  "The network integrates gene-pathway-module annotation relationships.",
  "Its edges are not experimentally validated physical interactions.",
  "Topological centrality does not establish biological causality.",
  "Network topology depends on the annotations and input gene signatures.",
  "Enrichment and network structure do not demonstrate pathway activation.",
  "",
  "OUTPUT FIGURES",
  "--------------",
  "01_full_network_overview.png",
  "02_top_betweenness_nodes.png",
  "03_bridge_gene_betweenness.png",
  "04_network_edge_composition.png",
  "05_bridge_gene_subnetwork.png",
  "06_community_size_distribution.png",
  "07_connected_component_sizes.png"
)

writeLines(
  report,
  file.path(report_dir, "network_visualization_interpretation.txt")
)

# ------------------------------------------------------------
# 20. AUDIT
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "input_gene_pathway_edges",
    "input_gene_module_edges",
    "input_pathway_module_edges",
    "unique_edges",
    "network_nodes",
    "network_edges",
    "bridge_genes",
    "connected_components",
    "communities",
    "top_nodes_exported"
  ),
  value = c(
    nrow(gene_pathway),
    nrow(gene_module),
    nrow(pathway_module),
    nrow(edges),
    vcount(graph),
    ecount(graph),
    length(bridge_names),
    length(component_result$csize),
    nrow(community_summary_39),
    nrow(top_30)
  )
)

write.csv(
  audit,
  file.path(audit_dir, "script39_audit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 21. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 39 COMPLETED\n")
cat("============================================================\n\n")

cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n")
cat("Connected components:", length(component_result$csize), "\n")
cat("Communities:", nrow(community_summary_39), "\n")
cat("Bridge genes:", length(bridge_names), "\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nFigures generated:\n")
cat(
  paste(
    list.files(figure_dir, pattern = "\\.png$"),
    collapse = "\n"
  ),
  "\n"
)

cat("\n============================================================\n")
