
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 38 - NETWORK TOPOLOGY AND FUNCTIONAL CONNECTIVITY
#
# Uses corrected Script 37 outputs.
# Network nodes: GENES, PATHWAYS, MODULES
# Network edges: GENE-PATHWAY, GENE-MODULE, PATHWAY-MODULE
#
# Analyses:
# - Degree centrality
# - Betweenness centrality
# - Harmonic centrality
# - Connected components
# - Articulation points
# - Louvain community detection
# - Bridge gene topology
#
# All network relationships represent annotations,
# NOT experimentally demonstrated molecular interactions.
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

PROJECT_ROOT <- normalizePath(
  "C:/Users/luank/Dropbox/PC/Desktop/PYTHON/TRANSCRIPTOME_X",
  winslash = "/",
  mustWork = FALSE
)

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_pathway_module_network",
  "tables"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "network_topology_analysis"
)

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
AUDIT_DIR <- file.path(OUTPUT_DIR, "audit")
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")

for (d in c(TABLE_DIR, AUDIT_DIR, REPORT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required_files <- c(
  "gene_pathway_edges.csv",
  "gene_module_edges.csv",
  "pathway_module_edges.csv",
  "gene_hub_ranking.csv",
  "multi_module_bridge_genes.csv"
)

for (f in required_files) {
  path <- file.path(INPUT_DIR, f)
  
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop(
    "Package 'igraph' is required. Install it with:\n",
    "install.packages('igraph')"
  )
}


# ============================================================
# 2. FUNCTIONS
# ============================================================

write_csv <- function(data, filename) {
  write.csv(
    data,
    filename,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}


# ============================================================
# 3. LOAD INPUTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 38 - NETWORK TOPOLOGY ANALYSIS\n")
cat("============================================================\n\n")

gene_pathway <- read.csv(
  file.path(INPUT_DIR, "gene_pathway_edges.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_module <- read.csv(
  file.path(INPUT_DIR, "gene_module_edges.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pathway_module <- read.csv(
  file.path(INPUT_DIR, "pathway_module_edges.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

hub_ranking <- read.csv(
  file.path(INPUT_DIR, "gene_hub_ranking.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bridge_genes <- read.csv(
  file.path(INPUT_DIR, "multi_module_bridge_genes.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. INPUT VALIDATION
# ============================================================

required_columns <- list(
  gene_pathway = c("gene_symbol", "pathway_id"),
  gene_module = c("gene_symbol", "module"),
  pathway_module = c("pathway_id", "module"),
  hub_ranking = c("gene_symbol", "hub_score", "hub_rank"),
  bridge_genes = c("gene_symbol", "module_count", "module_list")
)

datasets <- list(
  gene_pathway = gene_pathway,
  gene_module = gene_module,
  pathway_module = pathway_module,
  hub_ranking = hub_ranking,
  bridge_genes = bridge_genes
)

for (dataset_name in names(datasets)) {
  
  missing <- setdiff(
    required_columns[[dataset_name]],
    colnames(datasets[[dataset_name]])
  )
  
  if (length(missing) > 0) {
    stop(
      "Missing columns in ",
      dataset_name,
      ": ",
      paste(missing, collapse = ", ")
    )
  }
}

gene_pathway$gene_symbol <- clean_text(gene_pathway$gene_symbol)
gene_pathway$pathway_id <- clean_text(gene_pathway$pathway_id)

gene_module$gene_symbol <- clean_text(gene_module$gene_symbol)
gene_module$module <- clean_text(gene_module$module)

pathway_module$pathway_id <- clean_text(pathway_module$pathway_id)
pathway_module$module <- clean_text(pathway_module$module)

hub_ranking$gene_symbol <- clean_text(hub_ranking$gene_symbol)
hub_ranking$hub_score <- safe_numeric(hub_ranking$hub_score)
hub_ranking$hub_rank <- safe_numeric(hub_ranking$hub_rank)

bridge_genes$gene_symbol <- clean_text(bridge_genes$gene_symbol)
bridge_genes$module_count <- safe_numeric(bridge_genes$module_count)
bridge_genes$module_list <- clean_text(bridge_genes$module_list)

if ("hub_score" %in% colnames(bridge_genes)) {
  bridge_genes$hub_score <- safe_numeric(bridge_genes$hub_score)
}

if ("hub_rank" %in% colnames(bridge_genes)) {
  bridge_genes$hub_rank <- safe_numeric(bridge_genes$hub_rank)
}

if ("bridge_status" %in% colnames(bridge_genes)) {
  bridge_genes$bridge_status <- clean_text(bridge_genes$bridge_status)
}

gene_pathway <- unique(
  gene_pathway[
    nzchar(gene_pathway$gene_symbol) &
      nzchar(gene_pathway$pathway_id),
    ,
    drop = FALSE
  ]
)

gene_module <- unique(
  gene_module[
    nzchar(gene_module$gene_symbol) &
      nzchar(gene_module$module) &
      gene_module$module != "NO_DEFINED_MODULE",
    ,
    drop = FALSE
  ]
)

pathway_module <- unique(
  pathway_module[
    nzchar(pathway_module$pathway_id) &
      nzchar(pathway_module$module) &
      pathway_module$module != "NO_DEFINED_MODULE",
    ,
    drop = FALSE
  ]
)

hub_ranking <- hub_ranking[
  nzchar(hub_ranking$gene_symbol),
  ,
  drop = FALSE
]

hub_ranking <- hub_ranking[
  !duplicated(hub_ranking$gene_symbol),
  ,
  drop = FALSE
]

bridge_genes <- bridge_genes[
  nzchar(bridge_genes$gene_symbol),
  ,
  drop = FALSE
]

bridge_genes <- bridge_genes[
  !duplicated(bridge_genes$gene_symbol),
  ,
  drop = FALSE
]


# ============================================================
# 5. BUILD NETWORK EDGE TABLE
# ============================================================

edges_gp <- data.frame(
  from = paste0("G::", gene_pathway$gene_symbol),
  to = paste0("P::", gene_pathway$pathway_id),
  edge_type = "GENE_PATHWAY",
  weight = 1,
  stringsAsFactors = FALSE
)

edges_gm <- data.frame(
  from = paste0("G::", gene_module$gene_symbol),
  to = paste0("M::", gene_module$module),
  edge_type = "GENE_MODULE",
  weight = 1,
  stringsAsFactors = FALSE
)

edges_pm <- data.frame(
  from = paste0("P::", pathway_module$pathway_id),
  to = paste0("M::", pathway_module$module),
  edge_type = "PATHWAY_MODULE",
  weight = 1,
  stringsAsFactors = FALSE
)

edges <- unique(
  rbind(edges_gp, edges_gm, edges_pm)
)

if (nrow(edges) == 0) {
  stop("No network edges available.")
}


# ============================================================
# 6. BUILD NODE TABLE
# ============================================================

genes <- unique(sub("^G::", "", edges$from[grepl("^G::", edges$from)]))
genes <- unique(c(
  genes,
  sub("^G::", "", edges$to[grepl("^G::", edges$to)])
))

pathways <- unique(sub("^P::", "", edges$from[grepl("^P::", edges$from)]))
pathways <- unique(c(
  pathways,
  sub("^P::", "", edges$to[grepl("^P::", edges$to)])
))

modules <- unique(sub("^M::", "", edges$from[grepl("^M::", edges$from)]))
modules <- unique(c(
  modules,
  sub("^M::", "", edges$to[grepl("^M::", edges$to)])
))

vertices <- rbind(
  data.frame(
    name = paste0("G::", genes),
    node_type = "GENE",
    label = genes,
    stringsAsFactors = FALSE
  ),
  data.frame(
    name = paste0("P::", pathways),
    node_type = "PATHWAY",
    label = pathways,
    stringsAsFactors = FALSE
  ),
  data.frame(
    name = paste0("M::", modules),
    node_type = "MODULE",
    label = modules,
    stringsAsFactors = FALSE
  )
)

vertices <- unique(vertices)

cat("Genes:", length(genes), "\n")
cat("Pathways:", length(pathways), "\n")
cat("Modules:", length(modules), "\n")
cat("Edges:", nrow(edges), "\n\n")


# ============================================================
# 7. CREATE IGRAPH NETWORK
# ============================================================

network <- igraph::graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = vertices
)

network <- igraph::simplify(
  network,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = "first"
)

igraph::V(network)$node_type <- vertices$node_type[
  match(
    igraph::V(network)$name,
    vertices$name
  )
]

igraph::V(network)$label <- vertices$label[
  match(
    igraph::V(network)$name,
    vertices$name
  )
]

cat("Network created.\n")
cat("Vertices:", igraph::vcount(network), "\n")
cat("Edges:", igraph::ecount(network), "\n\n")


# ============================================================
# 8. CENTRALITY METRICS
# ============================================================

cat("Calculating degree centrality...\n")

degree_values <- igraph::degree(
  network,
  mode = "all"
)

cat("Calculating betweenness centrality...\n")

betweenness_values <- igraph::betweenness(
  network,
  directed = FALSE,
  normalized = TRUE
)

cat("Calculating harmonic centrality...\n")

harmonic_values <- igraph::harmonic_centrality(
  network,
  mode = "all",
  normalized = TRUE
)

components <- igraph::components(network)

articulation <- igraph::articulation_points(network)

articulation_names <- igraph::V(network)$name[
  articulation
]

community <- igraph::cluster_louvain(
  network,
  weights = igraph::E(network)$weight
)

community_membership <- igraph::membership(
  community
)


# ============================================================
# 9. NODE TOPOLOGY TABLE
# ============================================================

node_topology <- data.frame(
  node_id = igraph::V(network)$name,
  node_type = igraph::V(network)$node_type,
  label = igraph::V(network)$label,
  degree = as.numeric(degree_values),
  degree_centrality = as.numeric(
    degree_values / max(1, igraph::vcount(network) - 1)
  ),
  betweenness = as.numeric(betweenness_values),
  harmonic_centrality = as.numeric(harmonic_values),
  component_id = as.integer(components$membership),
  community_id = as.integer(community_membership),
  is_articulation_point = (
    igraph::V(network)$name %in% articulation_names
  ),
  stringsAsFactors = FALSE
)

node_topology <- node_topology[
  order(
    -node_topology$betweenness,
    -node_topology$degree
  ),
  ,
  drop = FALSE
]

rownames(node_topology) <- NULL


# ============================================================
# 10. GENE TOPOLOGY RANKING
# ============================================================

gene_topology <- node_topology[
  node_topology$node_type == "GENE",
  ,
  drop = FALSE
]

gene_topology$gene_symbol <- gene_topology$label

# Merge only hub-ranking fields.
# Module and bridge metadata are retained separately in bridge_genes.

hub_columns <- hub_ranking[
  ,
  c("gene_symbol", "hub_score", "hub_rank"),
  drop = FALSE
]

gene_topology <- merge(
  gene_topology,
  hub_columns,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

gene_topology <- gene_topology[
  order(
    -gene_topology$betweenness,
    -gene_topology$degree,
    -gene_topology$hub_score
  ),
  ,
  drop = FALSE
]

gene_topology$topology_rank <- seq_len(
  nrow(gene_topology)
)

rownames(gene_topology) <- NULL


# ============================================================
# 11. PATHWAY TOPOLOGY RANKING
# ============================================================

pathway_topology <- node_topology[
  node_topology$node_type == "PATHWAY",
  ,
  drop = FALSE
]

pathway_topology$pathway_id <- sub(
  "^P::",
  "",
  pathway_topology$node_id
)

pathway_topology <- pathway_topology[
  order(
    -pathway_topology$betweenness,
    -pathway_topology$degree
  ),
  ,
  drop = FALSE
]

pathway_topology$topology_rank <- seq_len(
  nrow(pathway_topology)
)


# ============================================================
# 12. MODULE TOPOLOGY
# ============================================================

module_topology <- node_topology[
  node_topology$node_type == "MODULE",
  ,
  drop = FALSE
]

module_topology$module <- module_topology$label

module_topology <- module_topology[
  order(
    -module_topology$betweenness,
    -module_topology$degree
  ),
  ,
  drop = FALSE
]


# ============================================================
# 13. BRIDGE GENE TOPOLOGY
# ============================================================

# Select topology-specific columns only.
# This prevents duplicate metadata columns and .x/.y suffixes.

topology_columns <- c(
  "gene_symbol",
  "node_id",
  "node_type",
  "label",
  "degree",
  "degree_centrality",
  "betweenness",
  "harmonic_centrality",
  "component_id",
  "community_id",
  "is_articulation_point",
  "topology_rank"
)

topology_columns <- intersect(
  topology_columns,
  colnames(gene_topology)
)

gene_topology_bridge <- gene_topology[
  ,
  topology_columns,
  drop = FALSE
]

bridge_topology <- merge(
  bridge_genes,
  gene_topology_bridge,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

# Add hub score/rank only if they are not already present.
missing_hub_columns <- setdiff(
  c("hub_score", "hub_rank"),
  colnames(bridge_topology)
)

if (length(missing_hub_columns) > 0) {
  
  hub_additional <- hub_ranking[
    ,
    c("gene_symbol", missing_hub_columns),
    drop = FALSE
  ]
  
  bridge_topology <- merge(
    bridge_topology,
    hub_additional,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )
}

# Ensure all columns required by the report exist.
required_bridge_output <- c(
  "gene_symbol",
  "module_count",
  "module_list",
  "degree",
  "betweenness",
  "harmonic_centrality",
  "community_id",
  "hub_score"
)

for (column_name in required_bridge_output) {
  
  if (!column_name %in% colnames(bridge_topology)) {
    bridge_topology[[column_name]] <- NA
  }
}

bridge_topology <- bridge_topology[
  order(
    -bridge_topology$betweenness,
    -bridge_topology$degree
  ),
  ,
  drop = FALSE
]

rownames(bridge_topology) <- NULL


# ============================================================
# 14. COMMUNITY ANALYSIS
# ============================================================

community_ids <- sort(
  unique(node_topology$community_id)
)

community_records <- lapply(
  community_ids,
  function(id) {
    
    d <- node_topology[
      node_topology$community_id == id,
      ,
      drop = FALSE
    ]
    
    data.frame(
      community_id = id,
      total_nodes = nrow(d),
      gene_nodes = sum(d$node_type == "GENE"),
      pathway_nodes = sum(d$node_type == "PATHWAY"),
      module_nodes = sum(d$node_type == "MODULE"),
      articulation_nodes = sum(d$is_articulation_point),
      mean_degree = round(mean(d$degree), 3),
      mean_betweenness = round(mean(d$betweenness), 6),
      stringsAsFactors = FALSE
    )
  }
)

community_summary <- do.call(
  rbind,
  community_records
)

community_summary <- community_summary[
  order(
    -community_summary$total_nodes
  ),
  ,
  drop = FALSE
]


# ============================================================
# 15. CONNECTED COMPONENT ANALYSIS
# ============================================================

component_ids <- sort(
  unique(node_topology$component_id)
)

component_records <- lapply(
  component_ids,
  function(id) {
    
    d <- node_topology[
      node_topology$component_id == id,
      ,
      drop = FALSE
    ]
    
    data.frame(
      component_id = id,
      total_nodes = nrow(d),
      gene_nodes = sum(d$node_type == "GENE"),
      pathway_nodes = sum(d$node_type == "PATHWAY"),
      module_nodes = sum(d$node_type == "MODULE"),
      stringsAsFactors = FALSE
    )
  }
)

component_summary <- do.call(
  rbind,
  component_records
)

component_summary <- component_summary[
  order(
    -component_summary$total_nodes
  ),
  ,
  drop = FALSE
]


# ============================================================
# 16. EDGE TYPE SUMMARY
# ============================================================

edge_summary <- as.data.frame(
  table(edges$edge_type),
  stringsAsFactors = FALSE
)

colnames(edge_summary) <- c(
  "edge_type",
  "edge_count"
)

edge_summary <- edge_summary[
  order(-edge_summary$edge_count),
  ,
  drop = FALSE
]


# ============================================================
# 17. NETWORK SUMMARY
# ============================================================

network_summary <- data.frame(
  metric = c(
    "total_nodes",
    "gene_nodes",
    "pathway_nodes",
    "module_nodes",
    "total_edges",
    "gene_pathway_edges",
    "gene_module_edges",
    "pathway_module_edges",
    "connected_components",
    "louvain_communities",
    "articulation_points",
    "multi_module_bridge_genes"
  ),
  value = c(
    igraph::vcount(network),
    sum(node_topology$node_type == "GENE"),
    sum(node_topology$node_type == "PATHWAY"),
    sum(node_topology$node_type == "MODULE"),
    igraph::ecount(network),
    sum(edges$edge_type == "GENE_PATHWAY"),
    sum(edges$edge_type == "GENE_MODULE"),
    sum(edges$edge_type == "PATHWAY_MODULE"),
    components$no,
    length(unique(community_membership)),
    length(articulation_names),
    nrow(bridge_topology)
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 18. EXPORT RESULTS
# ============================================================

write_csv(
  node_topology,
  file.path(TABLE_DIR, "node_topology_metrics.csv")
)

write_csv(
  gene_topology,
  file.path(TABLE_DIR, "gene_topology_ranking.csv")
)

write_csv(
  pathway_topology,
  file.path(TABLE_DIR, "pathway_topology_ranking.csv")
)

write_csv(
  module_topology,
  file.path(TABLE_DIR, "module_topology_ranking.csv")
)

write_csv(
  bridge_topology,
  file.path(TABLE_DIR, "bridge_gene_topology.csv")
)

write_csv(
  community_summary,
  file.path(TABLE_DIR, "network_community_summary.csv")
)

write_csv(
  component_summary,
  file.path(TABLE_DIR, "connected_component_summary.csv")
)

write_csv(
  edge_summary,
  file.path(TABLE_DIR, "network_edge_summary.csv")
)

write_csv(
  network_summary,
  file.path(TABLE_DIR, "network_topology_summary.csv")
)

write_csv(
  data.frame(
    node_id = names(community_membership),
    community_id = as.integer(community_membership),
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "node_community_assignments.csv")
)


# ============================================================
# 19. AUDIT
# ============================================================

audit <- data.frame(
  metric = c(
    "input_gene_pathway_edges",
    "input_gene_module_edges",
    "input_pathway_module_edges",
    "network_vertices",
    "network_edges",
    "gene_nodes",
    "pathway_nodes",
    "module_nodes",
    "connected_components",
    "louvain_communities",
    "articulation_points",
    "bridge_genes",
    "igraph_version"
  ),
  value = c(
    nrow(gene_pathway),
    nrow(gene_module),
    nrow(pathway_module),
    igraph::vcount(network),
    igraph::ecount(network),
    length(genes),
    length(pathways),
    length(modules),
    components$no,
    length(unique(community_membership)),
    length(articulation_names),
    nrow(bridge_topology),
    as.character(utils::packageVersion("igraph"))
  ),
  stringsAsFactors = FALSE
)

write_csv(
  audit,
  file.path(AUDIT_DIR, "script38_audit.csv")
)


# ============================================================
# 20. REPORT
# ============================================================

report_file <- file.path(
  REPORT_DIR,
  "network_topology_report.txt"
)

sink(report_file)

cat("============================================================\n")
cat("TRANSCRIPTOME-X\n")
cat("SCRIPT 38 - NETWORK TOPOLOGY ANALYSIS\n")
cat("============================================================\n\n")

cat("NETWORK SUMMARY\n")
print(network_summary, row.names = FALSE)

cat("\nEDGE TYPE SUMMARY\n")
print(edge_summary, row.names = FALSE)

cat("\nCOMMUNITY SUMMARY\n")
print(community_summary, row.names = FALSE)

cat("\nCONNECTED COMPONENTS\n")
print(component_summary, row.names = FALSE)

cat("\nTOP 20 GENES BY BETWEENNESS\n")
print(
  head(
    gene_topology[
      ,
      c(
        "gene_symbol",
        "degree",
        "betweenness",
        "harmonic_centrality",
        "community_id",
        "hub_score",
        "hub_rank"
      )
    ],
    20
  ),
  row.names = FALSE
)

cat("\nBRIDGE GENE TOPOLOGY\n")
print(
  bridge_topology[
    ,
    c(
      "gene_symbol",
      "module_count",
      "module_list",
      "degree",
      "betweenness",
      "harmonic_centrality",
      "community_id",
      "hub_score"
    )
  ],
  row.names = FALSE
)

cat("\nTOP 15 PATHWAYS BY BETWEENNESS\n")
print(
  head(
    pathway_topology[
      ,
      c(
        "pathway_id",
        "degree",
        "betweenness",
        "community_id"
      )
    ],
    15
  ),
  row.names = FALSE
)

cat("\nMODULE TOPOLOGY\n")
print(
  module_topology[
    ,
    c(
      "module",
      "degree",
      "betweenness",
      "harmonic_centrality",
      "community_id"
    )
  ],
  row.names = FALSE
)

cat("\nMETHODOLOGICAL NOTES\n")
cat(
  "- The network represents functional annotations and shared memberships.\n"
)
cat(
  "- Edges do not establish direct molecular interactions.\n"
)
cat(
  "- Betweenness measures shortest-path structural centrality.\n"
)
cat(
  "- Louvain communities represent topological groupings.\n"
)
cat(
  "- Community membership is not proof of biological pathway activation.\n"
)
cat(
  "- Hub scores from Script 37 and topology metrics are distinct measures.\n"
)

cat("\n============================================================\n")
cat("END OF SCRIPT 38\n")
cat("============================================================\n")

sink()


# ============================================================
# 21. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 38 COMPLETED\n")
cat("============================================================\n\n")

cat("Network nodes:", igraph::vcount(network), "\n")
cat("Network edges:", igraph::ecount(network), "\n")
cat("Connected components:", components$no, "\n")

cat(
  "Louvain communities:",
  length(unique(community_membership)),
  "\n"
)

cat("Articulation points:", length(articulation_names), "\n")

cat("\nTOP 10 GENES BY BETWEENNESS\n")

print(
  head(
    gene_topology[
      ,
      c(
        "gene_symbol",
        "degree",
        "betweenness",
        "community_id",
        "hub_score"
      )
    ],
    10
  ),
  row.names = FALSE
)

cat("\nBRIDGE GENES\n")

print(
  bridge_topology[
    ,
    c(
      "gene_symbol",
      "module_count",
      "degree",
      "betweenness",
      "community_id"
    )
  ],
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\n============================================================\n")
