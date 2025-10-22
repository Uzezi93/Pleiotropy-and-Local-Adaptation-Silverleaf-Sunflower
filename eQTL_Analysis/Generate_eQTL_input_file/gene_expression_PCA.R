# ---- Libraries ----
library(tidyverse)
library(cowplot)
library(ggfortify)

# ---- Paths ----
setwd("/Users/uzeziokinedo/Documents/Helianthus_argophyllus_project/Helianthus_argophyllus project/")
expr_path   <- "~/Downloads/expression_vst_matrix.txt"   # genes x samples, row.names=TRUE
coords_path <- "~/Downloads/arg_native_poplatlong.csv"   # columns: pop, long, lat, X (population label)

# ---- Load expression (genes x samples) ----
expr <- read.delim(expr_path, row.names = 1, check.names = FALSE)
stopifnot(ncol(expr) > 0, nrow(expr) > 0)

# ---- Load and tidy population info ----
coords_clean <- read.csv(coords_path, stringsAsFactors = FALSE) %>%
  distinct(pop, .keep_all = TRUE) %>%
  dplyr::select(pop, X, long, lat) %>%
  dplyr::rename(population = X)

# Define the accessions to re-label
gene_bank_samples <- c("Ames449", "Ames695", "ARG1805", "ARG1820", "ARG1834")

# Update population column for those samples
coords_clean <- coords_clean %>%
  mutate(population = ifelse(pop %in% gene_bank_samples, "Gene_Bank", population)) %>%
  filter(!population %in% c("Gene_Bank", "South"))

# Check that it worked
coords_clean %>%
  filter(pop %in% gene_bank_samples) %>%
  dplyr::select(pop, population)

#pop_sample <- coords_clean %>%
 # dplyr::select(pop, population) %>%
 # na.omit()

#names(pop_sample) <- c("sample", "population")
#write.table(pop_sample, file = "~/Downloads/keep.samples.txt", quote = FALSE, sep = "\t", row.names = F)

# sample metadata from expression column names
samples_all <- colnames(expr)
meta_tbl <- tibble(sample = samples_all) %>%
  mutate(prefix = sub("[.-].*$", "", sample)) %>%   # take part before first "." or "-"
  left_join(coords_clean, by = c("prefix" = "pop"))

if (anyNA(meta_tbl$population)) {
  cat("[WARN] Unmatched samples (dropped):",
      paste(meta_tbl$sample[is.na(meta_tbl$population)], collapse = ", "), "\n")
}

# ---- Keep *all* matched samples, including South ----
meta_filt <- meta_tbl %>%
  filter(!is.na(population))

keep_samples <- meta_filt$sample
expr <- expr[, intersect(colnames(expr), keep_samples), drop = FALSE]

# ensure alignment between expr and metadata
meta_filt <- meta_filt %>% dplyr::slice(match(colnames(expr), sample))
stopifnot(identical(colnames(expr), meta_filt$sample))
cat("[INFO] Kept", ncol(expr), "samples and", nrow(expr), "genes.\n")
table(meta_filt$population)

# ---- Prepare matrix for PCA (samples x genes) ----
expr_t <- t(as.matrix(expr))      # samples x genes
mode(expr_t) <- "numeric"

expr_t <- expr_t[rowSums(expr_t) > 0, , drop = FALSE]
keep_var <- apply(expr_t, 2, function(x) var(x, na.rm = TRUE) > 0)
keep_fin <- colSums(is.finite(expr_t)) == nrow(expr_t)
keep <- keep_var & keep_fin
cat("[INFO] Dropping", sum(!keep), "genes (zero-variance or NA/Inf).\n")
expr_t_use <- expr_t[, keep, drop = FALSE]

# align metadata to PCA samples
meta_for_plot <- meta_filt %>% dplyr::slice(match(rownames(expr_t_use), sample))
stopifnot(identical(rownames(expr_t_use), meta_for_plot$sample))

# ---- PCA (uncorrected) ----
pcat <- prcomp(expr_t_use, center = TRUE, scale. = TRUE)

custom_colors <- c("North" = "aquamarine4", "Coast" = "brown", "South" = "purple", "Gene_Bank" = "darkgrey")

uncorrected_PCA <- autoplot(
  pcat,
  data   = meta_for_plot,
  colour = "population",
  size   = 3
) +
  ggtitle("Uncorrected PCA") +
  scale_color_manual(values = custom_colors) +
  theme_light()

# ---- Scree plot ----
explained_variance <- (pcat$sdev^2) / sum(pcat$sdev^2)
p_scree <- tibble(
  PC = seq_along(explained_variance),
  Var = explained_variance
) %>%
  ggplot(aes(PC, Var)) +
  geom_col(fill = "#377eb8", color = "white", width = 0.8) +  # blue bars
  geom_line(color = "#e41a1c", size = 1) +                    # red line
  geom_point(color = "#e41a1c", size = 2) +                   # red points
  labs(
    title = "Scree plot",
    x = "Principal Component",
    y = "Variance Explained"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("expression_screeplot.tiff", plot = p_scree,
       device = "tiff", width = 7, height = 4, dpi = 300)

# ---- Regress out PC1 + PC2 ----
regress_out_pcs <- function(M, k = 2, center = TRUE, scale. = TRUE) {
  stopifnot(is.matrix(M), ncol(M) > 0, nrow(M) > 1)
  pca <- prcomp(M, center = center, scale. = scale.)
  k <- min(k, ncol(pca$x))
  X <- cbind(Intercept = 1, pca$x[, seq_len(k), drop = FALSE])
  XtX_inv <- solve(crossprod(X))
  B <- XtX_inv %*% crossprod(X, M)
  fitted <- X %*% B
  resid  <- M - fitted
  list(residuals = resid, pca = pca, betas = B)
}

res_out <- regress_out_pcs(expr_t_use, k = 2, center = TRUE, scale. = TRUE)
expr_resid <- res_out$residuals

# ---- PCA on residuals ----
pcat_corr <- prcomp(expr_resid, center = TRUE, scale. = TRUE)
corrected_PCA <- autoplot(
  pcat_corr,
  data   = meta_for_plot,
  colour = "population",
  size   = 3
) +
  ggtitle("Corrected PCA (PC1–PC2 regressed out)") +
  scale_color_manual(values = custom_colors) +
  theme_light()

expr_PCA <- plot_grid(uncorrected_PCA, corrected_PCA, nrow = 1, labels = c("A", "B"))
ggsave("expression_PCA.tiff", plot = expr_PCA, device = "tiff", width = 10, height = 5, dpi = 300)

# ---- Save matrices ----
write.table(expr_t_use, "expr_used_for_pca.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

write.table(expr_resid, "~/Downloads/expr_residual_matrix.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

write.table(t(expr_t_use), "expr_used_for_pca_gene_x_sample.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(t(expr_resid), "expr_residual_gene_x_sample.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("[INFO] Done. Wrote:\n",
    " - expression_screeplot.tiff\n",
    " - expression_PCA.tiff\n",
    " - expr_used_for_pca.txt\n",
    " - expr_residual_matrix.txt\n")

# ---- Save PC1 and PC2 covariates ----
pc_scores <- as.data.frame(pcat$x[, 1:2])
pc_scores <- tibble::rownames_to_column(pc_scores, var = "sample")

cov_mat <- t(pc_scores[, -1])
colnames(cov_mat) <- pc_scores$sample
rownames(cov_mat) <- c("PC1", "PC2")

write.table(cov_mat, file = "~/Downloads/covariates_pc12.txt", sep = "\t", quote = FALSE)



