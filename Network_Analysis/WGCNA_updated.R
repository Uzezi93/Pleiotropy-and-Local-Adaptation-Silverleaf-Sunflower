# -------------------------- Libraries --------------------------
suppressPackageStartupMessages({
  library(WGCNA)        
  library(tidyverse)    
  library(igraph)       
  library(scales)
  library(ggpubr)
})
allowWGCNAThreads()     

# -------------------------- Load expression --------------------
vst <- read.table("~/Downloads/expression_vst_matrix.txt",
                  header = TRUE, sep = "\t", check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)

# first column as rownames
first_col_name <- names(vst)[1]
if (!any(grepl("^[A-Za-z].*", first_col_name)) && ncol(vst) > 2) {
  # Heuristic fallback (rare)
  message("Using first column as gene IDs (heuristic).")
  rownames(vst) <- vst[[1]]
  vst[[1]] <- NULL
} else if (grepl("gene|Gene|Geneid|Gene_id", first_col_name)) {
  rownames(vst) <- vst[[1]]
  vst[[1]] <- NULL
}

# Ensure all entries are numeric and samples are rows, genes are columns
vst <- as.data.frame(vst)

# fix sample names 
names(vst) <- make.names(names(vst), unique = TRUE)

# transpose matrix and convert to numeric
datExpr <- t(as.matrix(vst))           # transpose: samples in rows
mode(datExpr) <- "numeric"

# Vector of sample names to drop
to_remove <- c("btm5.1", "arg11B.11", "btm7B.14",
               "Ames449", "Ames695", "ARG1805", "ARG1820", "ARG1834")

# Keep only rows that are not in to_remove
datExpr <- datExpr[!(rownames(datExpr) %in% to_remove), ]

# Check dimensions and row names after filtering
dim(datExpr)
head(rownames(datExpr))


# Basic QC: drop bad samples/genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# -------------------------- Soft-threshold ---------------------
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow = c(1,2)); cex1 <- 0.9
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     main="Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",
     type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")

# Choose a power (from your scree: 3–4 looked good; you chose 4)
picked_power <- 4

# -------------------------- Network construction ---------------
# Keep WGCNA's cor so nothing masks it
cor <- WGCNA::cor

netwk <- blockwiseModules(
  datExpr,
  power              = picked_power,
  networkType        = "signed",
  deepSplit          = 2,
  pamRespectsDendro  = FALSE,
  minModuleSize      = 30,
  maxBlockSize       = 4000,
  reassignThreshold  = 0,
  mergeCutHeight     = 0.25,
  saveTOMs           = TRUE,
  saveTOMFileBase    = "ER",
  numericLabels      = TRUE,
  verbose            = 3
)

# -------------------------- Quick diagnostics ------------------
plotDendroAndColors(
  netwk$dendrograms[[1]],
  labels2colors(netwk$colors[netwk$blockGenes[[1]]]),
  "Module colors (no merging)",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)


mergedColors <- labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[ netwk$blockGenes[[1]] ],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
table(mergedColors)

# -------------------------- Module eigengenes ------------------
ME_list <- moduleEigengenes(datExpr, colors = mergedColors, excludeGrey = TRUE, impute = TRUE)
MEs     <- orderMEs(ME_list$eigengenes)  # columns like "MEblue", "MEgreen", ...
write.table(MEs, "~/Downloads/module_eigengenes.tsv", sep="\t", quote=FALSE)


# -------------------------- Connectivity metrics --------------


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)   
  library(readr)
})

datExpr

# ---- 1) Connectivity (one value per gene/column) ----
connectivity <- softConnectivity(
  datExpr,                 # rows = samples, cols = genes
  corFnc     = "cor",
  weights    = NULL,
  type       = "signed",
  power      = 4,
  blockSize  = 4000,
  minNSamples= NULL,
  verbose    = 2, indent = 0
)

# Save with proper gene IDs (columns of datExpr)
connectivity_df <- tibble(
  gene         = colnames(datExpr),
  connectivity = as.numeric(connectivity)
)
write.table(connectivity_df, file = "~/Downloads/connectivity.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

# ---- 2) Load & normalize connectivity ----
connectivity2 <- read.table("~/Downloads/connectivity.txt", header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE) %>% as_tibble()

# Normalize gene IDs: remove "gene:" prefix and any ".<suffix>"
normalize_gene <- function(x) {
  x <- as.character(x)
  x <- sub("^gene:", "", x)
  x <- sub("\\..*$", "", x)
  x
}
connectivity2 <- connectivity2 %>%
  mutate(gene = normalize_gene(gene)) %>%
  distinct(gene, .keep_all = TRUE)

# ---- 3) Function to read gene lists robustly into a `gene` column ----
read_genes_any <- function(path, col = NULL) {
  df <- suppressWarnings(read.table(path, header = FALSE, sep = "\t",
                                    stringsAsFactors = FALSE, quote = "", comment.char = ""))
  if (!is.null(col)) {
    stopifnot(col <= ncol(df))
    genes <- df[[col]]
  } else {
    # auto-pick: prefer col with "gene:" content or longest string column
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
  tibble(gene = normalize_gene(genes)) %>% filter(!is.na(gene), gene != ".", gene != "")
}

# ---- 4) Load the sets (edit paths if needed) ----
lfmm_genes   <- read_genes_any("~/Downloads/lfmm_genes2.txt") %>% distinct()
pcadapt_genes<- read_genes_any("~/Downloads/pcadapt_genes2.txt") %>% distinct()
eQTL_genes   <- read_genes_any("~/Downloads/eQTL_genes2.txt")      %>% distinct()
eGene_genes  <- read_genes_any("~/Downloads/eGenes_genes2.txt")               %>% distinct()
control_genes <- read_genes_any("~/Downloads/control_genes2.txt") %>% distinct()
shared_loci <- read_genes_any("~/Downloads/common_outliers_genes2.txt") %>% distinct()

# ---- 5) Slice connectivity for each group ----
lfmm_connectivity    <- semi_join(connectivity2, lfmm_genes,   by = "gene")
pcadapt_connectivity <- semi_join(connectivity2, pcadapt_genes,by = "gene")
eQTL_connectivity    <- semi_join(connectivity2, eQTL_genes,   by = "gene")
eGene_connectivity   <- semi_join(connectivity2, eGene_genes,  by = "gene")
shared_connectivity <- semi_join(connectivity2, shared_loci,  by = "gene")

# Define Control as genes not in any of the outlier/eQTL/eGene sets:
out_union <- bind_rows(lfmm_genes, pcadapt_genes, eQTL_genes, eGene_genes) %>%
  distinct(gene)
control_connectivity <- anti_join(connectivity2, out_union, by = "gene")

# Save connectivity for each group
write.table(lfmm_connectivity,    "~/Downloads/lfmm_connectivity.txt",    sep="\t", quote=FALSE, row.names=FALSE)
write.table(pcadapt_connectivity, "~/Downloads/pcadapt_connectivity.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(eQTL_connectivity,    "~/Downloads/eQTL_connectivity.txt",    sep="\t", quote=FALSE, row.names=FALSE)
write.table(eGene_connectivity,   "~/Downloads/eGene_connectivity.txt",   sep="\t", quote=FALSE, row.names=FALSE)
write.table(control_connectivity, "~/Downloads/control_connectivity.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(shared_connectivity, "~/Downloads/shared_connectivity.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---- 6) Build one long table for plotting & tests ----
to_long <- function(df, label) df %>% transmute(Method = label, connectivity)
data_long <- bind_rows(
  to_long(eGene_connectivity,   "eGene"),
  to_long(eQTL_connectivity,    "eQTL"),
  to_long(lfmm_connectivity,    "LFMM"),
  to_long(pcadapt_connectivity, "PCAdapt"),
  to_long(control_connectivity, "Control")
) %>%
  mutate(connectivity = as.numeric(connectivity)) %>%
  filter(is.finite(connectivity))

# Wilcoxon tests vs Control (for reference)
con_control <- filter(data_long, Method == "Control")$connectivity
con_lfmm    <- filter(data_long, Method == "LFMM")$connectivity
con_pcad    <- filter(data_long, Method == "PCAdapt")$connectivity
print(wilcox.test(con_control, con_lfmm))
print(wilcox.test(con_control, con_pcad))
print(wilcox.test(con_pcad,    con_lfmm))

# Median line (Control)
m_ctrl <- median(con_control, na.rm = TRUE)

# ---- 7) Plot ----
p <- ggplot(data_long, aes(x = Method, y = connectivity, fill = Method)) +
  geom_violin(adjust = 2, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     ref.group = "Control") +
  scale_y_continuous(trans = "log10") +
  geom_hline(yintercept = m_ctrl, color = "black", alpha = 0.7) +
  labs(x = NULL, y = "log10(connectivity)", title = NULL) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = "connectivity.png", plot = p, width = 5, height = 4.5, dpi = 300)


# Check the enrichment of genes identified by both pcadapt and lfmm analysis in WGCNA network modules

## ---------- Helpers ----------
normalize_gene <- function(x){
  x <- as.character(x)
  x <- sub("^gene:", "", x)
  x <- sub("\\..*$", "", x)
  x
}

## ---------- Build modules data.frame aligned to datExpr ----------
# mergedColors must be the color names returned by labels2colors(netwk$colors)
modules1 <- tibble(
  gene   = normalize_gene(colnames(datExpr)),
  colors = as.character(mergedColors)
) %>%
  distinct(gene, .keep_all = TRUE)

# Drop "grey" module since they are unassigned modules
mods_keep <- dplyr::filter(modules1, colors != "grey")

## ---------- Define outlier sets (normalize IDs) ----------
lfmm_outliers    <- normalize_gene(lfmm_genes$gene)      |> unique()
pcadapt_outliers <- normalize_gene(pcadapt_genes$gene)   |> unique()

# Choose which set to test: union, intersection, etc.
outlier_genes <- union(lfmm_outliers, pcadapt_outliers)  # union
# outlier_genes <- intersect(lfmm_outliers, pcadapt_outliers)  # shared only

## ---------- Fisher’s exact test per module (enrichment) ----------
genes_u   <- mods_keep$gene
N         <- length(genes_u)
out_set   <- intersect(outlier_genes, genes_u)
tot_out   <- length(out_set)
stopifnot(tot_out > 0)

# count per module
tab_mod <- mods_keep %>%
  mutate(is_out = as.integer(gene %in% out_set)) %>%
  group_by(colors) %>%
  summarise(
    n_module      = n(),
    n_out_in_mod  = sum(is_out),
    .groups = "drop"
  ) %>%
  mutate(
    n_out_outside = tot_out - n_out_in_mod,
    n_non_in_mod  = n_module - n_out_in_mod,
    n_non_outside = (N - n_module) - n_out_outside
  )
write.table(fisher_results, "~/Downloads/module_fisher_union.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# run Fisher (one-sided 'greater' for enrichment), and 2-sided CI
fisher_results <- tab_mod %>%
  rowwise() %>%
  mutate(
    m = list(matrix(c(n_out_in_mod,  n_out_outside,
                      n_non_in_mod,  n_non_outside),
                    nrow = 2, byrow = TRUE)),
    ft_greater  = list(fisher.test(m, alternative = "greater")),
    ft_twosided = list(fisher.test(m, alternative = "two.sided")),
    odds_ratio  = as.numeric(ft_twosided$estimate),   # same estimate either way
    p_value     = ft_greater$p.value,                 # enrichment p-value
    ci_lower    = ft_twosided$conf.int[1],
    ci_upper    = ft_twosided$conf.int[2]
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(colors, n_module,
         n_out_in_mod, n_non_in_mod,
         odds_ratio, ci_lower, ci_upper, p_value, p_adj) %>%
  filter(p_adj < 0.01)

fisher_results_all <- tab_mod %>%
  rowwise() %>%
  mutate(
    m = list(matrix(c(n_out_in_mod,  n_out_outside,
                      n_non_in_mod,  n_non_outside),
                    nrow = 2, byrow = TRUE)),
    ft_greater  = list(fisher.test(m, alternative = "greater")),
    ft_twosided = list(fisher.test(m, alternative = "two.sided")),
    odds_ratio  = as.numeric(ft_twosided$estimate),   # same estimate either way
    p_value     = ft_greater$p.value,                 # enrichment p-value
    ci_lower    = ft_twosided$conf.int[1],
    ci_upper    = ft_twosided$conf.int[2]
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(colors, n_module,
                n_out_in_mod, n_non_in_mod,
                odds_ratio, ci_lower, ci_upper, p_value, p_adj)

write.table(fisher_results_all, "~/Downloads/module_fisher_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



## ---------- Permutation test (label randomization) ----------
# Keep module sizes fixed; shuffle which genes are called "outliers".
# This tests whether observed clustering of outliers into modules
# exceeds what I’d get by chance given module sizes.

permute_module_enrichment <- function(mods_df, outlier_genes, B = 10000, seed = 1){
  mods_df <- distinct(mods_df, gene, colors)
  genes   <- mods_df$gene
  N       <- length(genes)
  out_set <- intersect(unique(outlier_genes), genes)
  tot_out <- length(out_set)
  stopifnot(tot_out > 0)
  
  # membership matrix (genes x modules)
  mod_levels <- unique(mods_df$colors)
  M <- sapply(mod_levels, function(cl) as.integer(mods_df$colors == cl))
  colnames(M) <- mod_levels
  in_mod <- colSums(M)
  
  # observed indicator
  y_obs <- as.integer(genes %in% out_set)
  
  # Haldane–Anscombe adjusted OR (avoids Inf when any cell is 0)
  OR_fun <- function(y){
    a <- as.numeric(crossprod(y, M))      # outliers in module
    b <- tot_out - a                      # outliers outside
    c <- in_mod - a                       # non-outliers in module
    d <- (N - in_mod) - b                 # non-outliers outside
    ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  }
  
  # observed ORs
  OR_obs <- OR_fun(y_obs)
  a_obs  <- as.numeric(crossprod(y_obs, M))
  
  set.seed(seed)
  null_OR <- replicate(B, {
    y <- integer(N); y[sample.int(N, tot_out)] <- 1L
    OR_fun(y)
  })
  
  # empirical one-sided p: Pr(null OR >= observed OR)
  emp_p  <- (rowSums(sweep(null_OR, 1, OR_obs, `>=`)) + 1) / (B + 1)
  emp_fdr <- p.adjust(emp_p, method = "BH")
  
  tibble(
    colors        = mod_levels,
    n_module_size = in_mod,
    n_outliers    = a_obs,
    odds_ratio    = OR_obs,
    emp_p         = emp_p,
    emp_fdr       = emp_fdr,
    null_mean     = rowMeans(null_OR),
    null_sd       = apply(null_OR, 1, sd),
    null_q95      = apply(null_OR, 1, function(x) quantile(x, 0.95))
  ) %>% arrange(emp_p)
}

# Run permutation for the same outlier set you used above
res_perm <- permute_module_enrichment(mods_keep, outlier_genes, B = 10000, seed = 42)
write.table(res_perm, "~/Downloads/module_perm_union.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Quick comparison: Fisher vs permutation
res_combined <- fisher_results_all %>%
  dplyr::select(colors, fisher_or = odds_ratio, fisher_p = p_value, fisher_fdr = p_adj) %>%
  left_join(res_perm %>% dplyr::select(colors, perm_or = odds_ratio, perm_p = emp_p, perm_fdr = emp_fdr),
            by = "colors") %>%
  arrange(perm_p)

print(res_combined %>% head(15))

write.table(res_combined, "~/Downloads/module_perm_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# ===================== Significant modules → gene lists + unenriched controls =====================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(purrr) })

# 0) Safety checks
stopifnot(exists("mods_keep"), is.data.frame(mods_keep))
mods_keep <- mods_keep %>% distinct(gene, colors)

# 1) Pick significant modules (prefer permutation, else Fisher)
sig_cols <- character(0)

if (exists("res_perm") && is.data.frame(res_perm) && all(c("colors","emp_fdr") %in% names(res_perm))) {
  sig_cols <- res_perm %>%
    filter(is.finite(emp_fdr), emp_fdr < 0.05) %>%
    arrange(emp_fdr) %>% pull(colors) %>% unique()
} else if (exists("fisher_results") && is.data.frame(fisher_results) &&
           all(c("colors","p_adj") %in% names(fisher_results))) {
  sig_cols <- fisher_results %>%
    filter(is.finite(p_adj), p_adj < 0.05) %>%
    arrange(p_adj) %>% pull(colors) %>% unique()
}

cat("Significant modules: ", ifelse(length(sig_cols)==0, "none", paste(sig_cols, collapse = ", ")), "\n", sep = "")

# 2) Genes in significant modules
sig_gene_table <- tibble(gene = character(), colors = character())
if (length(sig_cols) > 0) {
  sig_gene_table <- mods_keep %>%
    filter(colors %in% sig_cols) %>%
    arrange(colors, gene)
}

# 3) Write gene lists (combined + per module)
if (nrow(sig_gene_table) > 0) {
  out_dir <- "~/Downloads/enriched_modules_gene_lists"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  readr::write_tsv(sig_gene_table, file.path(out_dir, "enriched_modules_genes.tsv"))
  
  sig_gene_table %>%
    group_by(colors) %>%
    group_walk(~ readr::write_lines(.x$gene,
                                    file.path(out_dir, paste0("genes_in_", unique(.x$colors), ".txt"))),
               .keep = TRUE)
  
  cat("Wrote enriched-module gene lists to: ", out_dir, "\n", sep = "")
} else {
  message("No significant modules to export.")
}

# 4) Build control pool from UNENRICHED modules and sample
set.seed(42)

unenriched_mods <- setdiff(unique(mods_keep$colors), sig_cols)
ctrl_pool <- mods_keep %>%
  filter(colors %in% unenriched_mods) %>%
  distinct(gene, .keep_all = TRUE)

n_enriched <- nrow(sig_gene_table)
if (n_enriched > 0) {
  if (n_enriched > nrow(ctrl_pool)) {
    warning("Requested control sample exceeds pool size; sampling with replacement.")
    ctrl_sample <- ctrl_pool %>%
      slice(sample.int(n(), size = n_enriched, replace = TRUE)) %>%
      mutate(sampled_with_replacement = TRUE)
  } else {
    ctrl_sample <- ctrl_pool %>%
      slice(sample.int(n(), size = n_enriched, replace = FALSE)) %>%
      mutate(sampled_with_replacement = FALSE)
  }
} else {
  ctrl_sample <- ctrl_pool %>% mutate(sampled_with_replacement = NA)
}

# per-module matched control: distribute enriched counts across unenriched modules ∝ size
permod_ctrl <- list()
if (n_enriched > 0) {
  permod_counts <- sig_gene_table %>% count(colors, name = "n_in_enriched")
  unenriched_sizes <- ctrl_pool %>% count(colors, name = "n_pool") %>% mutate(w = n_pool / sum(n_pool))
  
  if (nrow(unenriched_sizes) > 0) {
    for (i in seq_len(nrow(permod_counts))) {
      need_i <- permod_counts$n_in_enriched[i]
      draws <- floor(unenriched_sizes$w * need_i)
      remainder <- need_i - sum(draws)
      if (remainder > 0) {
        frac <- (unenriched_sizes$w * need_i) - draws
        add_ix <- order(frac, decreasing = TRUE)[seq_len(remainder)]
        draws[add_ix] <- draws[add_ix] + 1L
      }
      samp_i <- map2_dfr(unenriched_sizes$colors, draws, function(cl, k) {
        if (k <= 0) return(tibble(gene = character(), colors = character()))
        pool_cl <- ctrl_pool %>% filter(colors == cl)
        tibble(
          gene = sample(pool_cl$gene, size = k, replace = k > nrow(pool_cl)),
          colors = cl
        )
      }) %>% mutate(enriched_color_matched = permod_counts$colors[i])
      permod_ctrl[[i]] <- samp_i
    }
  }
}
ctrl_sample_permodule <- bind_rows(permod_ctrl)

# 5) Save control sets
out_dir_ctrl <- "~/Downloads/enriched_modules_controls"
dir.create(out_dir_ctrl, showWarnings = FALSE, recursive = TRUE)

readr::write_tsv(ctrl_pool,   file.path(out_dir_ctrl, "control_pool_unenriched_all.tsv"))          # full pool
readr::write_tsv(ctrl_sample, file.path(out_dir_ctrl, "control_unenriched_size_matched.tsv"))      # pooled size-matched
if (nrow(ctrl_sample_permodule) > 0) {
  readr::write_tsv(ctrl_sample_permodule,
                   file.path(out_dir_ctrl, "control_unenriched_permodule_matched.tsv"))            # per-module matched
}

cat("Unenriched control pool: ", nrow(ctrl_pool), " genes\n", sep = "")
cat("Size-matched pooled control: ", nrow(ctrl_sample), " genes (target = ", n_enriched, ")\n", sep = "")
if (nrow(ctrl_sample_permodule) > 0) {
  cat("Per-module matched control: ", nrow(ctrl_sample_permodule), " rows (sum over matched sets)\n", sep = "")
}

#-----------------------------------------------------------
# Get edge lists for cytoscape visualization

modules_of_interest <- c("blue","red","yellow", "pink", "green")  # example
keepGenes <- netwk$colors %in% match(modules_of_interest, labels2colors(1:20))

# Export network only for those
TOM <- TOMsimilarityFromExpr(datExpr[, keepGenes], power = picked_power, networkType = "signed")

exportNetworkToCytoscape(
  TOM,
  edgeFile = "CytoscapeInput-edges-selected.txt",
  nodeFile = "CytoscapeInput-nodes-selected.txt",
  weighted = TRUE,
  threshold = 0.05,
  nodeNames = colnames(datExpr)[keepGenes],
  nodeAttr = netwk$colors[keepGenes]
)


edge_list <- read.table("~/Documents/Helianthus_argophyllus_project/Helianthus_argophyllus project/CytoscapeInput-edges-selected.txt", skip = 1)
edge_list$V5 <-  NULL
edge_list$V6 <-  NULL
names(edge_list) <- c("From", "To", "Weight", "Direction")

edge_list2 <- edge_list %>%
  mutate(
    From = str_remove(From, "^gene:"),
    To   = str_remove(To, "^gene:")
  )

lfmm_genes <- read.table("~/Downloads/lfmm_genes2.txt")
pcadapt_genes <- read.table("~/Downloads/pcadapt_genes2.txt")
outlier_genes <- c(lfmm_genes$V1, pcadapt_genes$V1)

blue <- read.table("~/Downloads/enriched_modules_gene_lists/genes_in_blue.txt")
# Subset edges where From is a blue gene
blue_edge_list <- edge_list2 %>%
  filter(From %in% blue$V1)

# Add outlier status to blue_edge_list
blue_edge_list <- blue_edge_list %>%
  dplyr::mutate(
    From_outlier = ifelse(From %in% outlier_genes, "Outlier", "Non-outlier"),
    To_outlier   = ifelse(To %in% outlier_genes, "Outlier", "Non-outlier")
  )

# Optional: a combined flag if *either end* is an outlier
blue_edge_list <- blue_edge_list %>%
  dplyr::mutate(
    Any_outlier = ifelse(From_outlier == "Outlier" | To_outlier == "Outlier",
                         "Outlier_edge", "Non-outlier_edge")
  )

# Quick check
table(blue_edge_list$From_outlier)
table(blue_edge_list$Any_outlier)

blue_outlier_edges <- blue_edge_list %>%
  filter(From %in% outlier_genes | To %in% outlier_genes) %>%
  filter(Weight > 0.3)

write.table(blue_outlier_edges, file = "~/Downloads/blue_edge_list.tsv", sep = "\t", quote = FALSE, row.names = F)

#----------------------------------------------------------
pink <- read.table("~/Downloads/enriched_modules_gene_lists/genes_in_pink.txt")

pink_edge_list <- edge_list2 %>%
  filter(From %in% pink$V1)

# Add outlier status to blue_edge_list
pink_edge_list <- pink_edge_list %>%
  dplyr::mutate(
    From_outlier = ifelse(From %in% outlier_genes, "Outlier", "Non-outlier"),
    To_outlier   = ifelse(To %in% outlier_genes, "Outlier", "Non-outlier")
  )

# Optional: a combined flag if *either end* is an outlier
pink_edge_list <- pink_edge_list %>%
  dplyr::mutate(
    Any_outlier = ifelse(From_outlier == "Outlier" | To_outlier == "Outlier",
                         "Outlier_edge", "Non-outlier_edge")
  )

# Quick check
table(pink_edge_list$From_outlier)
table(pink_edge_list$Any_outlier)

pink_outlier_edges <- pink_edge_list %>%
  filter(From %in% outlier_genes | To %in% outlier_genes) %>%
  filter(Weight > 0.1)

write.table(pink_outlier_edges, file = "~/Downloads/pink_edge_list.tsv", sep = "\t", quote = FALSE, row.names = F)

#----------------------------------------------------------
red <- read.table("~/Downloads/enriched_modules_gene_lists/genes_in_red.txt")
red_edge_list <- edge_list2 %>%
  filter(From %in% red$V1)

# Add outlier status to blue_edge_list
red_edge_list <- red_edge_list %>%
  dplyr::mutate(
    From_outlier = ifelse(From %in% outlier_genes, "Outlier", "Non-outlier"),
    To_outlier   = ifelse(To %in% outlier_genes, "Outlier", "Non-outlier")
  )

# Optional: a combined flag if *either end* is an outlier
red_edge_list <- red_edge_list %>%
  dplyr::mutate(
    Any_outlier = ifelse(From_outlier == "Outlier" | To_outlier == "Outlier",
                         "Outlier_edge", "Non-outlier_edge")
  )

# Quick check
table(red_edge_list$From_outlier)
table(red_edge_list$Any_outlier)

red_outlier_edges <- red_edge_list %>%
  filter(From %in% outlier_genes | To %in% outlier_genes) %>%
  filter(Weight > 0.2)

write.table(red_outlier_edges, file = "~/Downloads/red_edge_list.tsv", sep = "\t", quote = FALSE, row.names = F)


yellow <- read.table("~/Downloads/enriched_modules_gene_lists/genes_in_yellow.txt")
yellow_edge_list <- edge_list2 %>%
  filter(From %in% yellow$V1)
yellow_edge_list$Attribute <- "yellow"
write.table(yellow_edge_list, file = "~/Downloads/yellow_edge_list.tsv", sep = "\t", quote = FALSE)

green <- read.table("~/Downloads/enriched_modules_gene_lists/genes_in_green.txt")
green_edge_list <- edge_list2 %>%
  filter(From %in% green$V1)
green_edge_list$Attribute <- "green"
write.table(green_edge_list, file = "~/Downloads/green_edge_list.tsv", sep = "\t", quote = FALSE)

