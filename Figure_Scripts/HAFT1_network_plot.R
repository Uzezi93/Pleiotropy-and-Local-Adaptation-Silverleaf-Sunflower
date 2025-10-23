#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(igraph); library(ggraph); library(ggplot2)
  library(fs);     library(tidyr);  library(forcats)
})

msg <- function(...) cat(sprintf(...), "\n")

# ----------------------- CONFIG -----------------------
bed_file <- "~/Downloads/genes.sorted.bed"

edge_files <- list(
  blue   = "~/Downloads/blue_edge_list.tsv",
  green  = "~/Downloads/green_edge_list.tsv",
  red    = "~/Downloads/red_edge_list.tsv",
  pink   = "~/Downloads/pink_edge_list.tsv",
  yellow = "~/Downloads/yellow_edge_list.tsv"
)

# Non-outlier color per module (plot)
module_colors <- c(
  blue   = "#4DBBD5",
  green  = "#00A087",
  red    = "#D64F4F",
  pink   = "#E07BB3",
  yellow = "#F4B400"
)

# HaFT1 search / window
haft1_pattern <- "HaFT1|FT1|Ha412HOChr06g"
haft1_chr    <- "Ha412HOChr06"
region_start <- 130000000L
region_end   <- 140000000L

min_edge_weight <- 0
max_nodes_label <- 120
out_dir <- "~/Downloads/out"
dir_create(out_dir, recurse = TRUE)

ensure_cols <- function(df, needed) {
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))
}

# If outlier_genes is absent, keep empty vector 
if (!exists("outlier_genes")) {
  msg("NOTE: 'outlier_genes' not found in environment; defaulting to empty vector.")
  outlier_genes <- character(0)
}

# ----------------------- READ BED & DEFINE WINDOW -----------------------
msg("Reading BED-like gene coordinates: %s", bed_file)
bed <- read.table(bed_file, header = FALSE, sep = "", check.names = FALSE,
                  quote = "", comment.char = "")
colnames(bed) <- c("chrom","start","end","gene","dot","strand")
bed$gene <- sub("^gene:", "", bed$gene)
bed <- bed %>%
  mutate(start = as.integer(start), end = as.integer(end),
         chrom = as.character(chrom), gene = as.character(gene))

msg("Locating HaFT1 by pattern: /%s/", haft1_pattern)
haft1_hits <- bed %>% filter(str_detect(gene, haft1_pattern))
if (nrow(haft1_hits) > 0) {
  msg("Found %d gene entries matching HaFT1 pattern.", nrow(haft1_hits))
  haft1_row <- haft1_hits %>%
    arrange(chrom, start) %>%
    { if ("Ha412HOChr06" %in% .$chrom) filter(., chrom == "Ha412HOChr06") else . } %>%
    slice(1)
  
  haft1_chr    <- haft1_row$chrom[1]
  haft1_start  <- haft1_row$start[1]
  haft1_end    <- haft1_row$end[1]
  flank_kb     <- 3000
  region_start <- max(0L, haft1_start - flank_kb*1000L)
  region_end   <- haft1_end + flank_kb*1000L
  msg("HaFT1 @ %s:%s-%s; window %s:%s-%s",
      haft1_chr, format(haft1_start, big.mark=","), format(haft1_end, big.mark=","),
      haft1_chr, format(region_start, big.mark=","), format(region_end, big.mark=","))
} else {
  msg("No gene matched HaFT1 pattern; using fixed window %s:%s-%s",
      haft1_chr, format(region_start, big.mark=","), format(region_end, big.mark=","))
}

region_genes <- bed %>%
  filter(chrom == haft1_chr, start >= region_start, end <= region_end) %>%
  distinct(gene, .keep_all = TRUE)

write_csv(region_genes, file.path(out_dir, "haft1_region_genes.csv"))
msg("Region genes: %d", nrow(region_genes))
stopifnot(nrow(region_genes) > 0)

region_gene_ids <- region_genes$gene

# Direct outliers limited to window (unique)
outlier_genes_window <- intersect(unique(outlier_genes), region_gene_ids)
msg("Outlier genes provided (unique): %d; in HaFT1 window: %d",
    length(unique(outlier_genes)), length(outlier_genes_window))

# ----------------------- PER-MODULE PROCESSOR -----------------------
process_module <- function(edge_path, module_name, module_color) {
  
  msg("\n=== Module: %s ===", module_name)
  msg("Reading edge list: %s", edge_path)
  
  edges <- read.table(edge_path, header = TRUE, sep = "", check.names = FALSE,
                      quote = "", comment.char = "")
  names(edges) <- gsub("\\s+", "", names(edges))
  ensure_cols(edges, c("From","To"))
  if (!"Weight" %in% names(edges)) edges$Weight <- 1.0
  
  # Keep edges that touch HaFT1-window genes
  region_edges <- edges %>%
    filter(From %in% region_gene_ids | To %in% region_gene_ids)
  if (min_edge_weight > 0)
    region_edges <- region_edges %>% filter(Weight >= min_edge_weight)
  
  write_csv(region_edges, file.path(out_dir, sprintf("haft1_region_edges_%s.csv", module_name)))
  msg("Region edges kept: %d", nrow(region_edges))
  stopifnot(nrow(region_edges) > 0)
  
  # Build graph (no extra isolated vertices)
  g <- graph_from_data_frame(region_edges, directed = FALSE)
  
  # ----- Node statuses for TABLE (two flags) -----
  # Direct outlier (your list)
  outlier_direct <- V(g)$name %in% outlier_genes_window
  
  # Partner-of-outlier (touches an outlier in this module), but not itself an outlier
  partners <- unique(c(
    region_edges$From[region_edges$To %in% outlier_genes_window],
    region_edges$To[ region_edges$From %in% outlier_genes_window]
  ))
  partner_of_outlier <- V(g)$name %in% setdiff(partners, outlier_genes_window)
  
  # Table 'outlier' column mirrors your original behavior (no propagation)
  node_tbl <- tibble(
    gene               = V(g)$name,
    in_region          = ifelse(V(g)$name %in% region_gene_ids, "Region", "Other"),
    is_HaFT1           = ifelse(str_detect(V(g)$name, haft1_pattern), "HaFT1", "Other"),
    degree             = degree(g),
    outlier            = ifelse(outlier_direct, "Outlier",
                                ifelse(V(g)$name %in% region_gene_ids, "Non-outlier", "Unknown")),
    outlier_direct     = outlier_direct,
    partner_of_outlier = partner_of_outlier
  )
  write_csv(node_tbl, file.path(out_dir, sprintf("haft1_region_nodes_%s.csv", module_name)))
  
  # ----- Labels: HaFT1 + its in-region neighbors + direct outliers -----
  haft1_ids <- V(g)$name[ str_detect(V(g)$name, haft1_pattern) ]
  if (length(haft1_ids) > 0) {
    nbr_ids <- unique(unlist(ego(g, order = 1, nodes = haft1_ids, mode = "all")))
    if (is.numeric(nbr_ids)) nbr_ids <- V(g)$name[nbr_ids]
  } else nbr_ids <- character(0)
  
  in_region_vec <- ifelse(V(g)$name %in% region_gene_ids, "Region", "Other")
  region_neighbors <- nbr_ids[in_region_vec[match(nbr_ids, V(g)$name)] == "Region"]
  
  label_set <- unique(c(haft1_ids, region_neighbors, V(g)$name[outlier_direct]))
  if (length(label_set) > max_nodes_label && max_nodes_label > 0) {
    ord  <- order(degree(g)[match(label_set, V(g)$name)], decreasing = TRUE)
    keep <- label_set[ord][seq_len(max_nodes_label)]
    V(g)$label <- ifelse(V(g)$name %in% keep, V(g)$name, "")
  } else {
    V(g)$label <- ifelse(V(g)$name %in% label_set, V(g)$name, "")
  }
  
  # ----- Plot: ONLY direct outliers get the outlier color; partners are NOT shown as outliers -----
  outlier_colors <- c(
    "Outlier"     = "brown4",          # direct outliers
    "Non-outlier" = module_color,      # module-matched base color
    "Unknown"     = "grey60"
  )
  plot_status <- ifelse(outlier_direct, "Outlier",
                        ifelse(V(g)$name %in% region_gene_ids, "Non-outlier", "Unknown"))
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, colour = "grey70", show.legend = FALSE) +
    geom_node_point(aes(size = pmax(1, log1p(degree(g)))),
                    colour = ifelse(in_region_vec == "Region", "grey85", NA),
                    stroke = 4, show.legend = FALSE) +
    geom_node_point(aes(color = plot_status,
                        shape = ifelse(str_detect(V(g)$name, haft1_pattern), "HaFT1", "Other"),
                        size  = pmax(1, log1p(degree(g)))),
                    stroke = 0.4) +
    geom_node_text(aes(label = V(g)$label),
                   size = 2.8, repel = TRUE, vjust = 1, family = "sans") +
    scale_color_manual(values = outlier_colors, name = "Gene status") +
    scale_shape_manual(values = c(HaFT1 = 17, Other = 16), name = "HaFT1") +
    scale_size(range = c(2, 8), guide = "none") +
    theme_void()
  
  png_file <- file.path(out_dir, sprintf("haft1_region_network_%s.png", module_name))
  pdf_file <- file.path(out_dir, sprintf("haft1_region_network_%s.pdf", module_name))
  ggsave(png_file, p, width = 10, height = 8, dpi = 300)
  ggsave(pdf_file, p, width = 10, height = 8)
  
  # Diagnostics
  tab_inwin <- table(factor(node_tbl$outlier[node_tbl$in_region == "Region"],
                            levels = c("Outlier","Non-outlier","Unknown")))
  msg("Saved plot & tables for %s.\nIn-window node labels:\n%s",
      module_name, paste(capture.output(print(tab_inwin)), collapse="\n"))
  
  invisible(list(graph=g, nodes=node_tbl, edges=region_edges, plot=p))
}

# ----------------------- RUN ALL MODULES -----------------------
results <- lapply(names(edge_files), function(m) {
  process_module(edge_files[[m]], module_name = m, module_color = module_colors[[m]])
})
names(results) <- names(edge_files)
msg("\nAll modules done. Outputs -> %s", out_dir)

