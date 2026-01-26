#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(cowplot)
})

allowWGCNAThreads()

# ============================================================
# CONFIG
# ============================================================
expr_path <- "~/Downloads/expression_vst_matrix.txt"   # expression matrix

out_dir <- "~/Downloads/outlier_enriched_module_centrality_size_matched"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# outlier sources (same as your pipeline)
lfmm_path    <- "~/Downloads/lfmm_genes2.txt"
pcadapt_path <- "~/Downloads/pcadapt_genes2.txt"

# Match your WGCNA settings
picked_power <- 4
network_type <- "signed"

# Same sample removals you used before
to_remove <- c("btm5.1", "arg11B.11", "btm7B.14",
               "Ames449", "Ames695", "ARG1805", "ARG1820", "ARG1834")

# Define enriched modules
fdr_cutoff_modules <- 0.05

# Size-matched centrality test settings (within each module)
B_perm <- 10000
seed   <- 42

# Combine top N module plots
n_panels_to_save <- 5

# Plot colors for some modules (others will fall back to the color name)
module_colors <- c(
  blue   = "#4DBBD5",
  green  = "#00A087",
  red    = "#D64F4F",
  pink   = "#E07BB3",
  yellow = "#F4B400"
)

# ============================================================
# Helpers
# ============================================================
normalize_gene <- function(x) {
  x <- as.character(x)
  x <- sub("^gene:", "", x)
  x <- sub("\\..*$", "", x)
  x
}

get_mod_color <- function(module_name) {
  if (!is.null(module_colors[[module_name]])) return(module_colors[[module_name]])
  return(module_name) # WGCNA color names usually work directly in ggplot
}

read_gene_list_lines <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x)]
  unique(normalize_gene(x))
}

# robust reader used earlier in your pipeline style (auto-picks column)
read_genes_any <- function(path, col = NULL) {
  df <- suppressWarnings(read.table(path, header = FALSE, sep = "\t",
                                    stringsAsFactors = FALSE, quote = "", comment.char = ""))
  if (!is.null(col)) {
    stopifnot(col <= ncol(df))
    genes <- df[[col]]
  } else {
    if (ncol(df) == 1) {
      genes <- df[[1]]
    } else {
      has_gene_prefix <- sapply(df, function(v) any(grepl("^gene:", v)))
      if (any(has_gene_prefix)) genes <- df[[ which(has_gene_prefix)[1] ]] else {
        lens <- sapply(df, function(v) mean(nchar(as.character(v)), na.rm = TRUE))
        genes <- df[[ which.max(lens) ]]
      }
    }
  }
  tibble(gene = normalize_gene(genes)) %>%
    filter(!is.na(gene), gene != ".", gene != "") %>%
    distinct(gene, .keep_all = TRUE)
}

p_to_stars <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}

# ============================================================
# 1) Load expression matrix and format for WGCNA
# ============================================================
message("Reading expression matrix: ", expr_path)
vst <- read.table(expr_path, header = TRUE, sep = "\t", check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)

gene_col <- names(vst)[1]
rownames(vst) <- normalize_gene(vst[[gene_col]])
vst[[gene_col]] <- NULL

# WGCNA expects: samples in rows, genes in columns
datExpr <- t(as.matrix(vst))
mode(datExpr) <- "numeric"

# Fix sample names for safety
rownames(datExpr) <- make.names(rownames(datExpr), unique = TRUE)

# Drop specified samples
datExpr <- datExpr[!(rownames(datExpr) %in% to_remove), , drop = FALSE]

# Basic QC
gsg <- goodSamplesGenes(datExpr, verbose = 2)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Normalize gene IDs again after QC
colnames(datExpr) <- normalize_gene(colnames(datExpr))
colnames(datExpr) <- make.unique(colnames(datExpr))

message("Expression dims (samples x genes): ", nrow(datExpr), " x ", ncol(datExpr))

# ============================================================
# 2) Run WGCNA modules
# ============================================================
cor <- WGCNA::cor

message("Running blockwiseModules...")
netwk <- blockwiseModules(
  datExpr,
  power              = picked_power,
  networkType        = network_type,
  deepSplit          = 2,
  pamRespectsDendro  = FALSE,
  minModuleSize      = 30,
  maxBlockSize       = 4000,
  reassignThreshold  = 0,
  mergeCutHeight     = 0.25,
  saveTOMs           = FALSE,
  numericLabels      = TRUE,
  verbose            = 3
)

mergedColors <- labels2colors(netwk$colors)

mods <- tibble(
  gene   = normalize_gene(colnames(datExpr)),
  colors = as.character(mergedColors)
) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(!is.na(colors), colors != "grey")

write_tsv(mods, file.path(out_dir, "gene_module_colors.tsv"))
message("Non-grey module-assigned genes: ", nrow(mods))
message("Modules: ", paste(sort(unique(mods$colors)), collapse = ", "))

# ============================================================
# 3) Load outliers (LFMM ∪ PCAdapt)
# ============================================================
lfmm_genes    <- read_genes_any(lfmm_path) %>% distinct()
pcadapt_genes <- read_genes_any(pcadapt_path) %>% distinct()

outlier_genes <- union(lfmm_genes$gene, pcadapt_genes$gene) %>% unique()
out_set <- intersect(outlier_genes, mods$gene)

stopifnot(length(out_set) > 0)
message("Outliers present in expression+modules: ", length(out_set))

write_lines(out_set, file.path(out_dir, "outliers_used_union_normalized.txt"))

# ============================================================
# 4) Identify outlier-enriched modules (Fisher + BH)
# ============================================================
all_genes <- mods$gene
N <- length(all_genes)
tot_out <- sum(all_genes %in% out_set)

tab_mod <- mods %>%
  mutate(is_out = gene %in% out_set) %>%
  group_by(colors) %>%
  summarise(
    n_module     = n(),
    n_out_in_mod = sum(is_out),
    .groups = "drop"
  ) %>%
  mutate(
    n_out_outside = tot_out - n_out_in_mod,
    n_non_in_mod  = n_module - n_out_in_mod,
    n_non_outside = (N - n_module) - n_out_outside
  )

fisher_tbl <- tab_mod %>%
  rowwise() %>%
  mutate(
    m = list(matrix(c(n_out_in_mod,  n_out_outside,
                      n_non_in_mod,  n_non_outside),
                    nrow = 2, byrow = TRUE)),
    ft = list(fisher.test(m, alternative = "greater")),
    odds_ratio = as.numeric(ft$estimate),
    p_value    = ft$p.value
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(fisher_tbl, file.path(out_dir, "module_outlier_enrichment_fisher.tsv"))

sig_cols <- fisher_tbl %>%
  filter(is.finite(p_adj), p_adj < fdr_cutoff_modules, n_out_in_mod > 0) %>%
  pull(colors) %>%
  unique()

if (length(sig_cols) == 0) {
  stop("No outlier-enriched modules at FDR < ", fdr_cutoff_modules,
       ". Consider relaxing fdr_cutoff_modules or verify outliers.")
}

message("Outlier-enriched modules: ", paste(sig_cols, collapse = ", "))

# ============================================================
# 5) Size-matched centrality test within each enriched module
#    + plot: dashed median (size-matched non-outliers) + perm p stars
#    + plot labels: Outliers / Non-outliers
#    + y shown as log10 values (units on axis are log10)
# ============================================================
test_module_centrality_size_matched <- function(module_color) {
  
  genes_mod <- mods %>% filter(colors == module_color) %>% pull(gene) %>% unique()
  
  # within-module soft connectivity (computed only among genes in this module)
  k <- softConnectivity(
    datExpr[, genes_mod, drop = FALSE],
    corFnc = "cor",
    type   = network_type,
    power  = picked_power,
    blockSize = 4000,
    verbose = 0
  )
  
  df <- tibble(
    gene = genes_mod,
    k_within = as.numeric(k),
    is_outlier = gene %in% out_set
  ) %>% filter(is.finite(k_within))
  
  n_out <- sum(df$is_outlier)
  n_non <- sum(!df$is_outlier)
  
  if (n_out < 2 || n_non < 2) {
    return(list(
      stats = tibble(
        colors = module_color,
        n_module = nrow(df),
        n_outliers = n_out,
        n_nonoutliers = n_non,
        median_out = NA_real_,
        median_non_size_matched = NA_real_,
        null_median_mean = NA_real_,
        effect_median_diff = NA_real_,
        perm_p = NA_real_
      ),
      plot = NULL
    ))
  }
  
  non_vals <- df$k_within[!df$is_outlier]
  med_out  <- median(df$k_within[df$is_outlier], na.rm = TRUE)
  
  # Null distribution: median of size-matched samples from NON-OUTLIERS
  set.seed(seed)
  replace_flag <- n_out > length(non_vals)
  null_meds <- replicate(B_perm, {
    samp <- sample(non_vals, size = n_out, replace = replace_flag)
    median(samp, na.rm = TRUE)
  })
  
  perm_p <- (sum(null_meds >= med_out) + 1) / (B_perm + 1)
  null_mean <- mean(null_meds, na.rm = TRUE)
  eff <- med_out - null_mean
  
  # Plot sample: ONE size-matched draw for visualization
  samp_plot <- sample(non_vals, size = n_out, replace = replace_flag)
  med_non_matched <- median(samp_plot, na.rm = TRUE)
  
  df_plot <- bind_rows(
    tibble(group = "Outliers",     k = df$k_within[df$is_outlier]),
    tibble(group = "Non-outliers", k = samp_plot)
  ) %>%
    mutate(group = factor(group, levels = c("Outliers", "Non-outliers")))
  
  # ---- PLOT ONLY: show log10 values explicitly on axis (unit-scale numbers)
  eps <- 1e-8
  df_plot <- df_plot %>% mutate(k_log10 = log10(k + eps))
  med_line <- log10(med_non_matched + eps)
  
  mod_col <- get_mod_color(module_color)
  stars <- p_to_stars(perm_p)
  
  y_star <- max(df_plot$k_log10, na.rm = TRUE) + 0.08
  
  p <- ggplot(df_plot, aes(x = group, y = k_log10)) +
    geom_violin(trim = FALSE, fill = mod_col, alpha = 0.6) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    # dashed line = median of SIZE-MATCHED non-outliers (log10 scale)
    geom_hline(yintercept = med_line, linetype = "dashed", linewidth = 0.6) +
    labs(
      x = NULL,
      y = "within module connectivity (Log10)",
      title = NULL,
      subtitle = paste0("Perm p=", signif(perm_p, 3),
                        " (", stars, ") | n_out=", n_out,
                        " | Outlier median=", signif(med_out, 3),
                        " | Non-outlier median=", signif(med_non_matched, 3))
    ) +
    annotate("text", x = 1.5, y = y_star, label = stars, size = 6) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10))
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0("centrality_size_matched_", module_color, ".png")),
    plot = p, width = 7, height = 4.2, dpi = 300
  )
  
  # Save values
  write_tsv(df, file.path(out_dir, paste0("centrality_values_full_", module_color, ".tsv")))
  write_tsv(tibble(null_median = null_meds),
            file.path(out_dir, paste0("centrality_null_medians_", module_color, ".tsv")))
  
  stats_row <- tibble(
    colors = module_color,
    n_module = nrow(df),
    n_outliers = n_out,
    n_nonoutliers = n_non,
    median_out = med_out,
    median_non_size_matched = med_non_matched,
    null_median_mean = null_mean,
    effect_median_diff = eff,
    perm_p = perm_p
  )
  
  list(stats = stats_row, plot = p)
}

# Run per module and collect plots + stats
res_list <- map(sig_cols, test_module_centrality_size_matched)
names(res_list) <- sig_cols

stats_tbl <- map_dfr(res_list, "stats") %>%
  mutate(perm_fdr = p.adjust(perm_p, method = "BH")) %>%
  arrange(perm_p)

write_tsv(stats_tbl, file.path(out_dir, "module_outlier_centrality_size_matched_stats.tsv"))

# ============================================================
# 6) Save combined multi-panel figure (A–E): 2 columns
#     - remove titles/subtitles ONLY for combined plot
# ============================================================
top_mods <- stats_tbl %>%
  filter(is.finite(perm_p)) %>%
  arrange(perm_p) %>%
  slice_head(n = n_panels_to_save) %>%
  pull(colors)

plot_map <- set_names(map(res_list, "plot"), sig_cols)

plots_to_combine <- plot_map[top_mods]
plots_to_combine <- plots_to_combine[!vapply(plots_to_combine, is.null, logical(1))]

# remove subtitle/title ONLY in combined plot
plots_to_combine <- lapply(plots_to_combine, function(pp) {
  pp + theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank()
  )
})

if (length(plots_to_combine) > 0) {
  combined <- cowplot::plot_grid(
    plotlist   = plots_to_combine,
    labels     = LETTERS[seq_along(plots_to_combine)],
    label_size = 14,
    ncol       = 2
  )
  
  nrows <- ceiling(length(plots_to_combine) / 2)
  
  ggsave(
    filename = file.path(out_dir, "plots", "centrality_size_matched_combined_AtoE.png"),
    plot = combined, width = 14, height = 4.2 * nrows, dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, "plots", "centrality_size_matched_combined_AtoE.pdf"),
    plot = combined, width = 14, height = 4.2 * nrows
  )
}

  
  message("Done.")
  message("Saved:")
  message(" - gene_module_colors.tsv")
  message(" - outliers_used_union_normalized.txt")
  message(" - module_outlier_enrichment_fisher.tsv")
  message(" - module_outlier_centrality_size_matched_stats.tsv")
  message(" - plots/centrality_size_matched_*.png")
  message(" - plots/centrality_size_matched_combined_AtoE.(png|pdf)")
  