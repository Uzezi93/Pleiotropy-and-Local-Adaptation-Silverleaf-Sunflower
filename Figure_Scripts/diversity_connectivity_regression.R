## ── Packages ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(broom); library(tidyr); library(stringr)
})

setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision")

## ── Load data ──────────────────────────────────────────────────────────────
connectivity <- read.table("connectivity.txt", header = TRUE, check.names = FALSE)
expr         <- read.table("expression_vst_stats.txt", header = TRUE, check.names = FALSE)
coast_pi     <- read_tsv("coast_pi_watterson_by_module.tsv", show_col_types = FALSE)  %>% dplyr::mutate(pop = "coast")
north_pi     <- read_tsv("north_pi_watterson_by_module.tsv", show_col_types = FALSE)  %>% dplyr::mutate(pop = "north")

## ── Build pooled gene×pop table ────────────────────────────────────────────
outlier_mods <- c("lfmm","pcadapt","shared")

pi_all <- dplyr::bind_rows(coast_pi, north_pi)

gene_outlier_flag <- pi_all %>%
  dplyr::mutate(is_out = module %in% outlier_mods) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    group_raw = ifelse(any(is_out, na.rm = TRUE), "outlier", "control"),
    .groups   = "drop"
  )

genes_pooled <- pi_all %>%
  dplyr::group_by(gene, pop) %>%
  dplyr::summarise(
    theta_pi = mean(pi, na.rm = TRUE),
    theta_w  = mean(watterson_theta, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::left_join(gene_outlier_flag, by = "gene") %>%
  dplyr::left_join(dplyr::select(connectivity, gene, connectivity), by = "gene") %>%
  dplyr::left_join(dplyr::select(expr, gene, meanExpr, exprVar), by = "gene") %>%
  dplyr::mutate(
    group        = factor(group_raw, levels = c("control","outlier")),
    pop          = factor(pop),
    connectivity = as.numeric(connectivity),
    expr_mean    = as.numeric(meanExpr),
    expr_var     = as.numeric(exprVar)
  ) %>%
  dplyr::select(gene, pop, group, theta_pi, theta_w, connectivity, expr_mean, expr_var)

## ── Standardize predictors & responses WITHIN GROUP (for comparability) ────
scale_within_group <- function(df) {
  df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      theta_pi_sc    = as.numeric(scale(theta_pi)),
      theta_w_sc     = as.numeric(scale(theta_w)),
      connectivity_z = as.numeric(scale(connectivity)),
      expr_mean_z    = as.numeric(scale(expr_mean)),
      expr_var_z     = as.numeric(scale(expr_var))
    ) %>%
    dplyr::ungroup()
}
genes_z <- scale_within_group(genes_pooled)

## ── Fit per-group models ────────────────
fit_one <- function(df, response) {
  fml <- stats::as.formula(paste(response, "~ connectivity_z + expr_mean_z + expr_var_z + pop"))
  fit <- stats::lm(fml, data = df, na.action = stats::na.exclude)
  list(coef = broom::tidy(fit, conf.int = TRUE), model = broom::glance(fit))
}

# Split by group
dat_ctrl <- dplyr::filter(genes_z, group == "control")
dat_out  <- dplyr::filter(genes_z, group == "outlier")

# Fit models
fit_ctrl_pi <- fit_one(dat_ctrl, "theta_pi_sc")
fit_ctrl_w  <- fit_one(dat_ctrl, "theta_w_sc")
fit_out_pi  <- fit_one(dat_out,  "theta_pi_sc")
fit_out_w   <- fit_one(dat_out,  "theta_w_sc")

## ── Build a wide summary table  ─────────────────────────
nice_num <- function(x) sprintf("%.2f", x)
cell_fmt <- function(est, lo, hi) paste0(nice_num(est), " (", nice_num(lo), ", ", nice_num(hi), ")")

pull_terms <- function(fit) {
  fit$coef %>%
    dplyr::filter(term %in% c("connectivity_z","expr_mean_z","expr_var_z")) %>%
    dplyr::transmute(
      term = dplyr::recode(term,
                           connectivity_z = "Connectivity",
                           expr_mean_z    = "Expression level",
                           expr_var_z     = "Expression variance"),
      est = estimate, lo = conf.low, hi = conf.high
    )
}

tbl_ctrl_pi <- pull_terms(fit_ctrl_pi) %>%
  dplyr::mutate(`Control θπ` = cell_fmt(est, lo, hi)) %>%
  dplyr::select(term, `Control θπ`)

tbl_ctrl_w  <- pull_terms(fit_ctrl_w) %>%
  dplyr::mutate(`Control θW` = cell_fmt(est, lo, hi)) %>%
  dplyr::select(term, `Control θW`)

tbl_out_pi  <- pull_terms(fit_out_pi) %>%
  dplyr::mutate(`Selection outliers θπ` = cell_fmt(est, lo, hi)) %>%
  dplyr::select(term, `Selection outliers θπ`)

tbl_out_w   <- pull_terms(fit_out_w) %>%
  dplyr::mutate(`Selection outliers θW` = cell_fmt(est, lo, hi)) %>%
  dplyr::select(term, `Selection outliers θW`)

# Merge rows
coef_block <- tbl_ctrl_pi %>%
  dplyr::full_join(tbl_ctrl_w,  by = "term") %>%
  dplyr::full_join(tbl_out_pi,  by = "term") %>%
  dplyr::full_join(tbl_out_w,   by = "term") %>%
  dplyr::mutate(term = factor(term, levels = c("Connectivity","Expression level","Expression variance"))) %>%
  dplyr::arrange(term)

# R² rows (one per response × group)
r2_rows <- tibble::tibble(
  term     = "R^2",
  `Control θπ`            = sprintf("%.2f", fit_ctrl_pi$model$r.squared),
  `Control θW`            = sprintf("%.2f", fit_ctrl_w$model$r.squared),
  `Selection outliers θπ` = sprintf("%.2f", fit_out_pi$model$r.squared),
  `Selection outliers θW` = sprintf("%.2f", fit_out_w$model$r.squared)
)

final_tbl <- dplyr::bind_rows(coef_block, r2_rows) %>%
  dplyr::rename(` ` = term) %>%   # single space to mimic blank header
  dplyr::select(` `, `Control θπ`, `Control θW`, `Selection outliers θπ`, `Selection outliers θW`)

## ── Save & print ───────────────────────────────────────────────────────────
write.table(final_tbl,
            file = "regression_per_group_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

print(final_tbl, n = Inf)
