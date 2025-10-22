# ---- libs ----
library(readr)
library(dplyr)
library(DESeq2)
library(matrixStats)

# ---------------- Load counts ----------------
expr <- read.delim(
  "~/Downloads/counts_all (1).txt",
  sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = "", stringsAsFactors = FALSE
)

stopifnot("Geneid" %in% names(expr))
stopifnot(ncol(expr) > 1)

# counts matrix (genes x samples)
cts <- as.matrix(expr[ , -1])
rownames(cts) <- expr$Geneid

# Make sure everything is numeric
if (!all(apply(cts, 2, is.numeric))) {
  cts <- apply(cts, 2, function(x) as.numeric(as.character(x)))
  rownames(cts) <- expr$Geneid
}

# ---------------- Build colData (one row per sample) ----------------
samples <- colnames(cts)
coldata <- data.frame(
  dummy = rep(1L, length(samples)),
  row.names = samples,
  check.names = FALSE
)

# sanity check
stopifnot(identical(rownames(coldata), colnames(cts)))

# ---------------- DESeq2 object ----------------
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ 1)

# ---------------- Filter low expression ----------------
# keep genes with >=10 counts in >=3 samples 
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("[INFO] Kept ", nrow(dds), " genes after low-expression filter.")

# ---------------- Normalization ----------------
# VST 
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
expr_norm <- assay(vsd)    # genes x samples (log-like scale)

write.table(
  cbind(Gene = rownames(expr_norm), as.data.frame(expr_norm, check.names = FALSE)),
  "~/Downloads/expression_vst_matrix.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---------------- Gene-level expression summary (from VST) 
# expr_norm is genes x samples (from vst)
stopifnot(is.matrix(expr_norm))

expr_stats_vst <- data.frame(
  gene        = rownames(expr_norm),
  meanExpr    = rowMeans(expr_norm, na.rm = TRUE),
  exprVar     = matrixStats::rowVars(expr_norm, na.rm = TRUE),
  exprSD      = sqrt(matrixStats::rowVars(expr_norm, na.rm = TRUE)),
  check.names = FALSE
)

# Save stats
write.table(
  expr_stats_vst,
  "~/Downloads/expression_vst_stats.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---------------- stats from size-factor–normalized counts --------
# (log1p for variance on counts scale)
norm_cts <- counts(dds, normalized = TRUE)
log_norm_cts <- log1p(norm_cts)

expr_stats_norm <- data.frame(
  gene        = rownames(log_norm_cts),
  meanExpr_ln = rowMeans(log_norm_cts, na.rm = TRUE),
  exprVar_ln  = matrixStats::rowVars(log_norm_cts, na.rm = TRUE),
  exprSD_ln   = sqrt(matrixStats::rowVars(log_norm_cts, na.rm = TRUE)),
  check.names = FALSE
)

write.table(
  expr_stats_norm,
  "~/Downloads/expression_normcount_stats.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

#----------------------------------------------------------
# uniquely mapped reads
library(tidyverse)

# Mapping stats from step 1
mapstats <- read_tsv("~/Downloads/star_mapstats.tsv", show_col_types = FALSE)

# Derive sample prefix to match your coords (everything before first "." or "-")
mapstats <- mapstats %>%
  mutate(prefix = sub("[.-].*$", "", sample))

# Load population info
coords <- read.csv("~/Downloads/arg_native_poplatlong.csv", stringsAsFactors = FALSE)

coords_clean <- coords %>%
  distinct(pop, .keep_all = TRUE) %>%
  transmute(prefix = pop, population = X)   # X is your North/Coast/South label

# Join
ms_anno <- mapstats %>%
  left_join(coords_clean, by = "prefix")

# Quick sanity check
ms_anno %>% filter(is.na(population)) %>% dplyr::select(sample) %>% distinct()

# If you want to compare only North vs Coast (exclude South):
ms_nc <- ms_anno %>% filter(population %in% c("North","Coast"))

write.table(ms_nc, file = "~/Downloads/star_mapping_stats.txt", sep = "\t", quote = F, row.names = F)

# Summary stats
ms_nc %>%
  group_by(population) %>%
  summarise(
    n = n(),
    mean_uniq = mean(uniq_map_pct, na.rm = TRUE),
    sd_uniq   = sd(uniq_map_pct, na.rm = TRUE),
    median_uniq = median(uniq_map_pct, na.rm = TRUE)
  )

# Statistical test (use Wilcoxon for robustness with small n)
wilcox.test(uniq_map_pct ~ population, data = ms_nc)

# Plot
ggplot(ms_nc, aes(population, uniq_map_pct, fill = population)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.8) +
  labs(x = NULL, y = "Uniquely mapped reads (%)",
       title = "STAR uniquely mapped reads by population") +
  theme_bw() +
  theme(legend.position = "none")

