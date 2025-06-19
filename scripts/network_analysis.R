#### STRING Network and Louvain Cluster Analysis for ba36_mild DEG ####

# --- Install and load required packages ---
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("tidygraph", quietly = TRUE)) install.packages("tidygraph")
if (!requireNamespace("clustermole", quietly = TRUE)) BiocManager::install("clustermole")

library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(stringr)
library(clustermole) # For cell type annotation

# --- 1. DEG Table Preparation ---

# Load DEG table (assume CSV with semicolon delimiter)
deg_data <- read.csv("ba36_mild.csv", sep = ";")

# Get unique gene list for STRING upload
unique_genes <- unique(deg_data$ensembl_gene_id)
write.table(unique_genes, "unique_gene_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- 2. STRING Network Import ---

# Load STRING interactions (TSV, columns: #node1, node2, combined_score, ...)
string_data <- read_tsv("string_interactions_short-dlpc.tsv", show_col_types = FALSE)

# Filter by confidence threshold and select relevant columns
interactions <- string_data %>%
  filter(combined_score > 0.4) %>%
  rename(gene1 = `#node1`, gene2 = node2) %>%
  select(gene1, gene2, combined_score)

# Build undirected igraph network with edge weights
full_network <- graph_from_data_frame(
  d = interactions,
  directed = FALSE
)
E(full_network)$weight <- interactions$combined_score

cat("Nodes:", vcount(full_network), "\n")
cat("Edges:", ecount(full_network), "\n")

# --- 3. Network Metrics (for reporting) ---
cat("Network density:", edge_density(full_network), "\n")
cat("Average node degree:", mean(degree(full_network)), "\n")
cat("Network diameter:", diameter(full_network), "\n")
cat("Average shortest path:", mean_distance(full_network), "\n")
cat("Clustering coefficient:", transitivity(full_network, type = "average"), "\n")
components <- components(full_network)
cat("Connected components:", components$no, "\n")
cat("Degree assortativity:", assortativity_degree(full_network), "\n")

# --- 4. Annotate DEGs with stage information ---
# Mutate gene attributes based on differential contrasts
gene_attributes <- deg_data %>%
  mutate(
    mild_sig     = mild_vs_healthy.y,
    moderate_sig = moderate_vs_healthy.y,
    severe_sig   = severe_vs_healthy.y,
    logFC_mild   = if_else(mild_sig     != 0, mild_vs_healthy.x, NA_real_),
    logFC_mod    = if_else(moderate_sig != 0, moderate_vs_healthy.x, NA_real_),
    logFC_sev    = if_else(severe_sig   != 0, severe_vs_healthy.x, NA_real_),
    n_sig = (mild_sig != 0) + (moderate_sig != 0) + (severe_sig != 0),
    logFC_sig_mean = case_when(
      n_sig > 0 ~ (coalesce(logFC_mild, 0) + coalesce(logFC_mod, 0) + coalesce(logFC_sev, 0)) / n_sig,
      TRUE ~ NA_real_
    )
  )

# --- 5. Build subnetwork for genes of interest ---
genes_of_interest <- gene_attributes$hgnc_symbol
common_genes <- intersect(V(full_network)$name, genes_of_interest)
sub_network <- induced_subgraph(full_network, vids = common_genes)

# Attach attributes to nodes
match_idx <- match(V(sub_network)$name, gene_attributes$hgnc_symbol)
V(sub_network)$mild_sig       <- gene_attributes$mild_sig[match_idx]
V(sub_network)$moderate_sig   <- gene_attributes$moderate_sig[match_idx]
V(sub_network)$severe_sig     <- gene_attributes$severe_sig[match_idx]
V(sub_network)$logFC_mild     <- gene_attributes$logFC_mild[match_idx]
V(sub_network)$logFC_mod      <- gene_attributes$logFC_mod[match_idx]
V(sub_network)$logFC_sev      <- gene_attributes$logFC_sev[match_idx]
V(sub_network)$n_sig          <- gene_attributes$n_sig[match_idx]
V(sub_network)$logFC_sig_mean <- gene_attributes$logFC_sig_mean[match_idx]

# --- 6. Louvain Clustering ---
clusters <- cluster_louvain(sub_network)
mem <- membership(clusters)
V(sub_network)$cluster <- mem

cat("Communities before filtering:", length(clusters), "\n")
print(sizes(clusters))

# Filter to clusters with at least 5 nodes
cluster_sizes <- sizes(clusters)
significant_clusters <- as.numeric(names(cluster_sizes[cluster_sizes >= 5]))
significant_nodes <- V(sub_network)[cluster %in% significant_clusters]
sub_network <- induced_subgraph(sub_network, vids = significant_nodes)
clusters <- cluster_louvain(sub_network)
mem <- membership(clusters)
V(sub_network)$cluster <- mem

cat("Communities after filtering:", length(clusters), "\n")
print(sizes(clusters))

# --- 7. Save annotated cluster info ---
cluster_df <- data.frame(
  gene           = names(mem),
  cluster        = as.integer(mem),
  mild_sig       = V(sub_network)$mild_sig,
  moderate_sig   = V(sub_network)$moderate_sig,
  severe_sig     = V(sub_network)$severe_sig,
  logFC_mild     = V(sub_network)$logFC_mild,
  logFC_mod      = V(sub_network)$logFC_mod,
  logFC_sev      = V(sub_network)$logFC_sev,
  n_sig          = V(sub_network)$n_sig,
  logFC_sig_mean = V(sub_network)$logFC_sig_mean,
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file = "clusters_ba36.csv", row.names = FALSE)

# --- 8. Identify local and global network hubs ---
gene_degree <- degree(sub_network, mode = "all")
top_hubs <- sort(gene_degree, decreasing = TRUE)[1:100]
global_hub_genes <- names(top_hubs)

# Local hubs (top 5 by degree in each cluster, degree > 2)
local_hubs <- lapply(unique(V(sub_network)$cluster), function(cl) {
  subg <- induced_subgraph(sub_network, vids = V(sub_network)[cluster == cl])
  degs <- degree(subg)
  top_degs <- head(sort(degs, decreasing = TRUE), 5)
  genes <- names(top_degs)
  data.frame(
    cluster         = rep(cl, length(genes)),
    gene            = genes,
    degree          = as.integer(top_degs),
    logFC_sig_mean  = V(subg)$logFC_sig_mean[ match(genes, V(subg)$name) ],
    n_sig           = V(subg)$n_sig[          match(genes, V(subg)$name) ],
    mild_sig        = V(subg)$mild_sig[       match(genes, V(subg)$name) ],
    moderate_sig    = V(subg)$moderate_sig[   match(genes, V(subg)$name) ],
    severe_sig      = V(subg)$severe_sig[     match(genes, V(subg)$name) ],
    stringsAsFactors = FALSE
  )
})
local_hubs_df <- bind_rows(local_hubs) %>% filter(degree > 2)

# --- 9. Visualization: Whole network with hubs ---
tg <- as_tbl_graph(sub_network) %>%
  activate(nodes) %>%
  mutate(
    hub_type = case_when(
      name %in% global_hub_genes   ~ "Global hub",
      name %in% local_hubs_df$gene ~ "Local hub",
      TRUE                         ~ "None"
    ),
    sig_stages = str_trim(paste0(
      ifelse(mild_sig     != 0, "mild ", ""),
      ifelse(moderate_sig != 0, "mod ", ""),
      ifelse(severe_sig   != 0, "sev ", "")
    )),
    sig_stages = str_replace(sig_stages, ", $", "")
  )
set.seed(42)
ggraph(tg, layout = "fr") +
  geom_edge_link(color = "grey80", alpha = 0.5) +
  geom_node_point(aes(
    color = factor(cluster),
    size  = abs(logFC_sig_mean),
    shape = hub_type
  ), alpha = 0.9) +
  geom_node_text(aes(
    label = paste0(name, "\n", sig_stages),
    fontface = ifelse(hub_type != "None", "bold", "plain")
  ), repel = TRUE, size = 2.5, lineheight = 0.9, max.overlaps = 50) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  scale_size_continuous(range = c(3, 8), name = "|mean logFC|") +
  scale_shape_manual(
    values = c("Global hub" = 17, "Local hub" = 15, "None" = 16),
    name = "Hub type"
  ) +
  theme_void() +
  ggtitle("STRING Subgraph: Louvain Clusters and Significant Stages") +
  theme(legend.position = "right")

# --- 10. (Optional) Cell type annotation using clustermole ---

markers <- clustermole_markers(species = "hs")
gene2ct <- markers %>%
  filter(organ == "Brain") %>%
  group_by(gene) %>%
  summarise(
    cell_type_list = paste(sort(unique(celltype)), collapse = ";"),
    .groups        = "drop"
  ) %>%
  mutate(cell_type_main = sub(";.*", "", cell_type_list))

cluster_df <- cluster_df %>%
  left_join(gene2ct, by = "gene") %>%
  mutate(cell_type_main = replace_na(cell_type_main, "Unknown"))

# (Optional) Summarize cell type composition in each cluster
cluster_ct <- cluster_df %>%
  count(cluster, cell_type_main, sort = TRUE) %>%
  slice_max(n, by = cluster, with_ties = FALSE)
print(cluster_ct)

# --- 11. Focused visualization: clusters with mild or moderate significance ---
clusters_to_show <- unique(
  V(sub_network)$cluster[
    V(sub_network)$mild_sig     != 0 |
    V(sub_network)$moderate_sig != 0
  ]
)
sub_network_filt <- induced_subgraph(
  sub_network,
  vids = V(sub_network)[cluster %in% clusters_to_show]
)
tg_filt <- as_tbl_graph(sub_network_filt) %>%
  activate(nodes) %>%
  mutate(
    hub_type = case_when(
      name %in% global_hub_genes   ~ "Global hub",
      name %in% local_hubs_df$gene ~ "Local hub",
      TRUE                         ~ "None"
    ),
    sig_stages = str_trim(paste0(
      ifelse(mild_sig     != 0, "mild, ", ""),
      ifelse(moderate_sig != 0, "mod, " , "")
    )),
    sig_stages = str_replace(sig_stages, ", $", "")
  )
layout <- create_layout(tg_filt, layout = "fr")
set.seed(42)
ggraph(layout) +
  geom_edge_link(color = "grey80", alpha = 0.5) +
  geom_node_point(aes(
    x     = x, y = y,
    color = factor(cluster),
    size  = abs(logFC_sig_mean),
    shape = hub_type
  ), alpha = 0.9) +
  geom_node_text(
    data = filter(layout, hub_type != "None"),
    aes(x = x, y = y, label = paste0(name, "\n", sig_stages)),
    fontface    = "bold", repel = TRUE, size = 3,
    lineheight  = 0.9, max.overlaps = 50
  ) +
  geom_node_text(
    data = filter(layout, hub_type == "None" & (mild_sig != 0 | moderate_sig != 0)),
    aes(x = x, y = y, label = paste0(name, "\n", sig_stages)),
    fontface    = "plain", repel = TRUE, size = 3,
    lineheight  = 0.9, max.overlaps = 50
  ) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  scale_size_continuous(range = c(3, 8), name = "|mean logFC|") +
  scale_shape_manual(
    values = c("Global hub" = 17, "Local hub" = 15, "None" = 16),
    name = "Hub type"
  ) +
  theme_void() +
  ggtitle("STRING: Labels for Hubs and Nodes with mild/moderate significance") +
  theme(legend.position = "right")

# --- 12. Export all cluster results ---
write.csv(cluster_df, file = "clusters_ba36.csv", row.names = FALSE)
