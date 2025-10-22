#!/usr/bin/env Rscript

# =========================================================
# Visualize HaFT1-region genes on a coexpression network
# with node colors indicating Outlier vs Non-outlier genes
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(fs)
  library(tidyr)
  library(forcats)
})

# -----------------------
# CONFIG (edit these)
# -----------------------
edge_file <- "~/Downloads/blue_edge_list.tsv"   
bed_file  <- "~/Downloads/genes.sorted.bed"            

# Try to detect HaFT1 by name; fallback to a window if not found
haft1_pattern <- "HaFT1|FT1|Ha412HOChr06g"   

# HaFT1 window
haft1_chr   <- "Ha412HOChr06"
region_start <- 130000000L
region_end   <- 140000000L
msg("Using fixed HaFT1 window %s:%s-%s",
    haft1_chr, format(region_start, big.mark=","), format(region_end, big.mark=","))


# Graph filtering
min_edge_weight <- 0        # set >0 to prune weak edges (e.g., 0.05)
max_nodes_label <- 120      # cap labels to avoid clutter

# Output folder
out_dir <- "~/Downloads/out"

# -----------------------
# Functions
# -----------------------
msg <- function(...) cat(sprintf(...), "\n")

ensure_cols <- function(df, needed) {
  missing <- setdiff(needed, names(df))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse=", "))
}

# Treat anything not "Non-outlier" (case-insensitive) as Outlier
norm_outlier_flag <- function(x) {
  y <- ifelse(is.na(x), "Unknown", as.character(x))
  ifelse(grepl("^\\s*non[- ]?outlier\\s*$", y, ignore.case = TRUE), "Non-outlier",
         ifelse(y == "Unknown", "Unknown", "Outlier"))
}

# -----------------------
# Read data
# -----------------------
msg("Reading edge list: %s", edge_file)
edges <- read.table(edge_file, header = TRUE, sep = "", check.names = FALSE, quote = "", comment.char = "")

# Standardize column names
names(edges) <- gsub("\\s+", "", names(edges))
if (!all(c("From", "To") %in% names(edges))) {
  stop("Edge list must contain 'From' and 'To' columns.")
}
if (!"Weight" %in% names(edges)) edges$Weight <- 1.0

msg("Reading BED-like gene coordinates: %s", bed_file)
bed <- read.table(bed_file, header = FALSE, sep = "", check.names = FALSE, quote = "", comment.char = "")
colnames(bed) <- c("chrom", "start", "end", "gene", "dot", "strand")
bed$gene <- sub("^gene:", "", bed$gene)

bed <- bed %>%
  mutate(start = as.integer(start), end = as.integer(end),
         chrom = as.character(chrom), gene = as.character(gene))

# -----------------------
# Locate HaFT1 and window
# -----------------------
msg("Locating HaFT1 by pattern: /%s/", haft1_pattern)
haft1_hits <- bed %>% filter(str_detect(gene, haft1_pattern))

if (nrow(haft1_hits) > 0) {
  msg("Found %d gene entries matching HaFT1 pattern.", nrow(haft1_hits))
  haft1_row <- haft1_hits %>%
    arrange(chrom, start) %>%
    { if ("Ha412HOChr06" %in% .$chrom) filter(., chrom == "Ha412HOChr06") else . } %>%
    slice(1)
  
  haft1_chr   <- haft1_row$chrom[1]
  haft1_start <- haft1_row$start[1]
  haft1_end   <- haft1_row$end[1]
  
  flank_kb <- 3000  # 3 Mb each side
  region_start <- max(0, haft1_start - flank_kb*1000)
  region_end   <- haft1_end + flank_kb*1000
  msg("HaFT1 @ %s:%s-%s; window %s:%s-%s",
      haft1_chr, format(haft1_start, big.mark=","), format(haft1_end, big.mark=","),
      haft1_chr, format(region_start, big.mark=","), format(region_end, big.mark=","))
} else {
  msg("No gene matched HaFT1 pattern; using window.")
  haft1_chr   <- haft1_chr_fallback
  region_start <- max(0, haft1_window_center - haft1_window_halfkb*1000)
  region_end   <- haft1_window_center + haft1_window_halfkb*1000
  msg("window %s:%s-%s",
      haft1_chr, format(region_start, big.mark=","), format(region_end, big.mark=","))
}

# -----------------------
# Genes in region
# -----------------------
region_genes <- bed %>%
  filter(chrom == haft1_chr,
         start >= region_start, end <= region_end) %>%
  distinct(gene, .keep_all = TRUE)

if (!dir_exists(out_dir)) dir_create(out_dir)
write_csv(region_genes, file.path(out_dir, "haft1_region_genes.csv"))

msg("Region genes: %d found.", nrow(region_genes))
if (!nrow(region_genes)) stop("No genes found in the specified region.")

region_gene_ids <- region_genes$gene

# -----------------------
# Edges involving region genes
# -----------------------
region_edges <- edges %>%
  filter(From %in% region_gene_ids | To %in% region_gene_ids)

if (min_edge_weight > 0) {
  region_edges <- region_edges %>% filter(Weight >= min_edge_weight)
}

write_csv(region_edges, file.path(out_dir, "haft1_region_edges.csv"))
msg("Region edges: %d kept (min_edge_weight=%s).", nrow(region_edges), min_edge_weight)
if (!nrow(region_edges)) stop("No edges involving region genes after filtering.")

# -----------------------
# Build graph
# -----------------------
g <- graph_from_data_frame(region_edges, directed = FALSE)

# -----------------------
# NODE-LEVEL OUTLIER STATUS
#   Derive from From_outlier / To_outlier if available
# -----------------------
node_outlier <- NULL
if (all(c("From_outlier","To_outlier") %in% names(region_edges))) {
  node_outlier <- bind_rows(
    region_edges %>% select(Gene = From, Status = From_outlier),
    region_edges %>% select(Gene = To,   Status = To_outlier)
  ) %>%
    mutate(Status = norm_outlier_flag(Status)) %>%
    group_by(Gene) %>%
    summarise(Gene_outlier = ifelse(any(Status == "Outlier"), "Outlier",
                                    ifelse(all(Status == "Non-outlier"), "Non-outlier", "Unknown")),
              .groups = "drop")
} else {
  # If no per-endpoint columns, fall back to Unknown
  node_outlier <- tibble(Gene = V(g)$name, Gene_outlier = "Unknown")
}

# Attach node attributes
V(g)$in_region  <- ifelse(V(g)$name %in% region_gene_ids, "Region", "Other")
V(g)$is_HaFT1   <- ifelse(str_detect(V(g)$name, haft1_pattern), "HaFT1", "Other")
V(g)$degree     <- degree(g)
V(g)$outlier    <- node_outlier$Gene_outlier[match(V(g)$name, node_outlier$Gene)]
V(g)$outlier    <- ifelse(is.na(V(g)$outlier), "Unknown", V(g)$outlier)

# Save node table
node_tbl <- tibble(
  gene       = V(g)$name,
  in_region  = V(g)$in_region,
  is_HaFT1   = V(g)$is_HaFT1,
  degree     = V(g)$degree,
  outlier    = V(g)$outlier
)
write_csv(node_tbl, file.path(out_dir, "haft1_region_nodes.csv"))

# -----------------------
# Label control — label HaFT1 nodes + their in-region neighbors
# -----------------------

# 1) Identify HaFT1 nodes
haft1_ids <- V(g)$name[V(g)$is_HaFT1 == "HaFT1"]

# 2) Collect 1-hop neighbors of ALL HaFT1 nodes
if (length(haft1_ids) > 0) {
  # ego() returns a list of vertex sequences; unlist to names
  nbr_ids <- unique(unlist(ego(g, order = 1, nodes = haft1_ids, mode = "all")))
  # Convert vertex indices to names if needed
  if (is.numeric(nbr_ids)) nbr_ids <- V(g)$name[nbr_ids]
} else {
  nbr_ids <- character(0)
}

# 3) Keep only neighbors that are physically inside the HaFT1 region
in_region_flag <- V(g)$in_region[match(nbr_ids, V(g)$name)] == "Region"
region_neighbors <- nbr_ids[!is.na(in_region_flag) & in_region_flag]

# 4) Final label set = (all HaFT1 nodes) ∪ (HaFT1 neighbors within region)
label_set <- unique(c(haft1_ids, region_neighbors))

# 5) To avoid clutter (uses max_nodes_label if you set it > 0)
if (length(label_set) > max_nodes_label && max_nodes_label > 0) {
  ord <- order(V(g)$degree[match(label_set, V(g)$name)], decreasing = TRUE)
  keep <- label_set[ord][seq_len(max_nodes_label)]
  V(g)$label <- ifelse(V(g)$name %in% keep, V(g)$name, "")
} else {
  V(g)$label <- ifelse(V(g)$name %in% label_set, V(g)$name, "")
}


# -----------------------
# Plot: color by Outlier status, shape by HaFT1, halo for region
# -----------------------
set.seed(42)

# Color map for outlier status
outlier_colors <- c(
  "Outlier" = "brown3",      # red
  "Non-outlier" = "#4DBBD5",  # blue
  "Unknown" = "grey60"
)

p <- ggraph(g, layout = "fr") +
  # edges
  geom_edge_link(alpha = 0.25, colour = "grey70", show.legend = FALSE) +
  # region halo (drawn first so it appears behind)
  geom_node_point(data = as.data.frame(igraph::as_data_frame(g, what = "vertices")) %>%
                    mutate(x = NA, y = NA),
                  inherit.aes = FALSE) +
  # base nodes with halo for region genes
  geom_node_point(aes(size = pmax(1, log1p(degree))),
                  colour = ifelse(V(g)$in_region == "Region", "grey85", NA),
                  stroke = 4, show.legend = FALSE) +
  # colored nodes (outlier status)
  geom_node_point(aes(color = outlier,
                      shape = is_HaFT1,
                      size  = pmax(1, log1p(degree))),
                  stroke = 0.4) +
  # labels
  geom_node_text(aes(label = label),
                 size = 2.8, repel = TRUE, vjust = 1, family = "sans") +
  scale_color_manual(values = outlier_colors, name = "Gene status") +
  scale_shape_manual(values = c(HaFT1 = 17, Other = 16), name = "HaFT1") +
  scale_size(range = c(2, 8), guide = "none") +
  theme_void() 

if (!dir_exists(out_dir)) dir_create(out_dir)
png_file <- file.path(out_dir, "haft1_region_network.png")
pdf_file <- file.path(out_dir, "haft1_region_network.pdf")
ggsave(png_file, p, width = 10, height = 8, dpi = 300)
ggsave(pdf_file, p, width = 10, height = 8)

msg("Saved plot to:\n - %s\n - %s", png_file, pdf_file)
msg("Saved tables to:\n - %s\n - %s\n - %s",
    file.path(out_dir,"haft1_region_genes.csv"),
    file.path(out_dir,"haft1_region_edges.csv"),
    file.path(out_dir,"haft1_region_nodes.csv"))

msg("Top 10 region genes by degree:")
print(
  node_tbl %>%
    arrange(desc(degree)) %>%
    head(10)
)
