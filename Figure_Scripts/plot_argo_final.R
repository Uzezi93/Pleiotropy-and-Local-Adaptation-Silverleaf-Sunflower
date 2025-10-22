#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggpubr)
  library(cowplot)
})

setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision/")

#----------------------------------------------------------------------------------------
# Proportion of eQTLs and eGenes in selection outliers
## ---------- helper: read unique gene IDs (col 4) ----------
read_gene_ids <- function(path) {
  tb <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, quote = "")
  if (ncol(tb) < 4) stop("Expected at least 4 columns in: ", path)
  unique(tb[[4]])
}

## ---------- load sets as gene IDs (V4) ----------
lfmm_ids    <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/lfmm_genes.coords.bed")
pcadapt_ids <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/pcadapt_genes.coords.bed")
eQTL_ids    <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/eQTL_genes.coords.bed")
eGene_ids   <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/eGene_genes.coords.bed")
All_ids     <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/control_genes.coords.bed")
shared_ids  <- read_gene_ids("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/shared_genes.coords.bed")

## ---------- proportions & plot ----------
prop_bars_clean <- function(lfmm, pcadapt, All, target, ylab = "Percent") {
  lfmm    <- unique(as.character(lfmm))
  pcadapt <- unique(as.character(pcadapt))
  All     <- unique(as.character(All))
  target  <- unique(as.character(target))
  sets <- list(All = All, LFMM = lfmm, PCAdapt = pcadapt)
  
  df <- lapply(names(sets), function(g) {
    S <- sets[[g]]; n <- length(S); x <- sum(S %in% target)
    tibble(Group = g, n = n, x = x, Proportion = x / n)
  }) %>% bind_rows() %>%
    mutate(Group = factor(Group, levels = c("All", "LFMM", "PCAdapt")))
  
  fisher_L <- with(df, fisher.test(matrix(c(
    x[Group=="LFMM"],    n[Group=="LFMM"]    - x[Group=="LFMM"],
    x[Group=="All"],     n[Group=="All"]     - x[Group=="All"]
  ), nrow = 2, byrow = TRUE)))
  fisher_P <- with(df, fisher.test(matrix(c(
    x[Group=="PCAdapt"], n[Group=="PCAdapt"] - x[Group=="PCAdapt"],
    x[Group=="All"],     n[Group=="All"]     - x[Group=="All"]
  ), nrow = 2, byrow = TRUE)))
  labs_p <- paste0("LFMM p=", signif(fisher_L$p.value, 2),
                   " | PCAdapt p=", signif(fisher_P$p.value, 2))
  
  ggplot(df, aes(Group, Proportion, fill = Group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(x, "/", n)), vjust = -0.4, size = 3.3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       limits = c(0, max(df$Proportion) * 1.20)) +
    scale_fill_manual(values = c(All = "grey70", LFMM = "brown", PCAdapt = "aquamarine4")) +
    labs(y = ylab, x = NULL, subtitle = labs_p) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

p1 <- prop_bars_clean(lfmm_ids, pcadapt_ids, All_ids, eQTL_ids,  "eQTL genes (%)")
p2 <- prop_bars_clean(lfmm_ids, pcadapt_ids, All_ids, eGene_ids, "eGene genes (%)")

grid_both <- cowplot::plot_grid(p1, p2, labels = c("A","B"), label_size = 18, ncol = 2, align = "hv")
ggplot2::ggsave("CONTROL_LFMM_PCADAPT_eQTL_eGene_grid.tiff", grid_both, width = 11, height = 4.2, dpi = 300)

# ---------------- palette (display labels) ----------------
pal <- c(
  "Control" = "darkgrey",
  "eQTL"    = "brown",
  "eGene"   = "brown",
  "LFMM"    = "aquamarine4",
  "PCAdapt" = "aquamarine4",
  "Shared"  = "steelblue"
)
mod_order_disp <- c("Control","LFMM","PCAdapt","Shared","eGene","eQTL")

module_to_disp <- function(x) {
  dplyr::recode(tolower(x),
                "control" = "Control", "lfmm" = "LFMM", "pcadapt" = "PCAdapt",
                "egene" = "eGene", "egenes" = "eGene", "eqtl" = "eQTL", "shared" = "Shared",
                .default = stringr::str_to_title(x)
  )
}

# ---------------- bootstrap control (for plotting distributions) ----------------
bootstrap_control_match <- function(df, value_col, control_name="control", seed=123){
  set.seed(seed)
  df_ctl <- df %>% dplyr::filter(module == control_name) %>% tidyr::drop_na({{value_col}})
  df_oth <- df %>% dplyr::filter(module != control_name) %>% tidyr::drop_na({{value_col}})
  if (nrow(df_ctl) == 0 || nrow(df_oth) == 0) return(df %>% tidyr::drop_na({{value_col}}))
  
  n_per_mod <- df_oth %>% dplyr::count(module, name = "n_mod")
  ctl_pool  <- df_ctl %>% dplyr::select(dplyr::all_of(setdiff(names(df_ctl), "module")))
  template_cols <- names(df)
  
  ctl_boot <- n_per_mod %>%
    dplyr::group_by(module) %>%
    dplyr::group_modify(function(.x, .y){
      n_m <- .x$n_mod[1]; ctl_pool %>% dplyr::slice_sample(n = n_m, replace = TRUE)
    }) %>% dplyr::ungroup() %>%
    dplyr::mutate(module = control_name) %>%
    dplyr::select(dplyr::all_of(template_cols))
  
  dplyr::bind_rows(df_oth, ctl_boot)
}

# ---------------- stats (raw Wilcoxon + bootstrap-balanced) ----------------
wilcox_against_control <- function(df, value_col, control_name="control", reps=1000, seed=123){
  set.seed(seed)
  df <- df %>% dplyr::mutate(module_disp = module_to_disp(module)) %>% tidyr::drop_na({{value_col}})
  present <- intersect(mod_order_disp, unique(df$module_disp))
  tests   <- setdiff(present, "Control")
  if (length(tests) == 0 || !("Control" %in% present)) return(tibble())
  
  ctl_vals <- df %>% dplyr::filter(module_disp == "Control") %>% dplyr::pull({{value_col}})
  
  purrr::map_dfr(tests, function(tg){
    x <- df %>% dplyr::filter(module_disp == tg) %>% dplyr::pull({{value_col}})
    p_raw <- tryCatch(wilcox.test(ctl_vals, x, exact = FALSE)$p.value, error = function(e) NA_real_)
    if (length(ctl_vals) > 0 && length(x) > 0) {
      ps <- replicate(reps, {
        xb <- sample(ctl_vals, size = length(x), replace = TRUE)
        suppressWarnings(wilcox.test(xb, x, exact = FALSE)$p.value)
      })
      tibble(group = tg, n_test = length(x), p_raw = p_raw,
             p_boot_median = median(ps, na.rm = TRUE),
             p_boot_min    = min(ps, na.rm = TRUE),
             p_boot_max    = max(ps, na.rm = TRUE))
    } else {
      tibble(group = tg, n_test = length(x), p_raw = NA_real_,
             p_boot_median = NA_real_, p_boot_min = NA_real_, p_boot_max = NA_real_)
    }
  })
}

# ===================== π / θW (your existing style) =====================
coast <- read_tsv("coast_pi_watterson_by_module.tsv", show_col_types = FALSE)
north <- read_tsv("north_pi_watterson_by_module.tsv", show_col_types = FALSE)

mod_order <- c("control","eGene","eQTL","lfmm","pcadapt","shared")

coast_pi <- coast %>%
  dplyr::mutate(module = factor(module, levels = mod_order)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta)
north_pi <- north %>%
  dplyr::mutate(module = factor(module, levels = mod_order)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta)

coast_pi_boot  <- coast_pi %>% bootstrap_control_match(pi)
coast_wat_boot <- coast_pi %>% bootstrap_control_match(watterson_theta)
north_pi_boot  <- north_pi  %>% bootstrap_control_match(pi)
north_wat_boot <- north_pi  %>% bootstrap_control_match(watterson_theta)

plot_violin_box <- function(df, value_col, title, ylab_expr, outfile = NULL){
  dfp <- df %>% dplyr::mutate(module_disp = module_to_disp(module)) %>% tidyr::drop_na({{value_col}})
  mods_present <- intersect(mod_order_disp, unique(dfp$module_disp))
  dfp <- dfp %>% dplyr::mutate(module_disp = factor(module_disp, levels = mods_present))
  
  ctl_med <- dfp %>%
    dplyr::filter(module_disp == "Control") %>%
    dplyr::summarise(med = median({{value_col}}, na.rm = TRUE)) %>% dplyr::pull(med)
  
  ptab <- wilcox_against_control(df, {{value_col}}, reps = 1000) %>%
    dplyr::mutate(star = dplyr::case_when(
      is.na(p_raw) ~ "", p_raw < 0.001 ~ "***", p_raw < 0.01 ~ "**", p_raw < 0.05 ~ "*", TRUE ~ ""
    )) %>% dplyr::filter(star != "")
  
  y_pos <- dfp %>% dplyr::group_by(module_disp) %>%
    dplyr::summarise(y = stats::quantile({{value_col}}, 0.98, na.rm = TRUE) * 1.05, .groups = "drop")
  
  stars_df <- ptab %>% dplyr::transmute(module_disp = factor(group, levels = levels(dfp$module_disp)),
                                        label = star) %>% dplyr::left_join(y_pos, by = "module_disp")
  
  p <- ggplot(dfp, aes(x = module_disp, y = {{value_col}}, fill = module_disp)) +
    geom_violin(trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    { if (!is.na(ctl_med)) geom_hline(yintercept = ctl_med, linetype = "dashed") } +
    { if (nrow(stars_df) > 0) geom_text(data = stars_df, aes(x = module_disp, y = y, label = label),
                                        inherit.aes = FALSE, fontface = "bold", size = 5) } +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = pal[mods_present], drop = FALSE, guide = "none") +
    labs(x = NULL, y = ylab_expr, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1))
  
  if (!is.null(outfile)) ggsave(outfile, p, width = 9, height = 5.5, dpi = 300)
  invisible(list(plot = p, pvals = ptab))
}

ylab_theta_pi <- "pi"
ylab_theta_w  <- "watterson theta"

p_coast_pi  <- plot_violin_box(coast_pi_boot,  pi,              "Coast", ylab_theta_pi, "coast_pi_violin_box_boot.png")
p_coast_w   <- plot_violin_box(coast_wat_boot, watterson_theta, "Coast", ylab_theta_w,  "coast_watterson_violin_box_boot.png")
p_north_pi  <- plot_violin_box(north_pi_boot,  pi,              "North", ylab_theta_pi, "north_pi_violin_box_boot.png")
p_north_w   <- plot_violin_box(north_wat_boot, watterson_theta, "North", ylab_theta_w,  "north_watterson_violin_box_boot.png")

combined <- cowplot::plot_grid(
  p_coast_pi$plot, p_north_pi$plot,
  p_coast_w$plot,  p_north_w$plot,
  labels = c("A","B","C","D"), ncol = 2, align = "hv"
)
ggsave("pi_watterson_coast_north_cowplot.png", combined, width = 12, height = 10, dpi = 300)

cat("\n[Wilcoxon vs Control — COAST π]\n");  print(wilcox_against_control(coast_pi,  pi, reps = 1000))
cat("\n[Wilcoxon vs Control — COAST θw]\n"); print(wilcox_against_control(coast_pi,  watterson_theta, reps = 1000))
cat("\n[Wilcoxon vs Control — NORTH π]\n");  print(wilcox_against_control(north_pi,  pi, reps = 1000))
cat("\n[Wilcoxon vs Control — NORTH θw]\n"); print(wilcox_against_control(north_pi,  watterson_theta, reps = 1000))

#-----------------------------------------------------------------------------------------
# ====================== FST / dXY with separate plotters ======================
fst <- read_tsv("fst_by_module_from_bed.tsv", show_col_types = FALSE)
dxy <- read_tsv("dxy_by_module_from_bed.tsv", show_col_types = FALSE)

fst_df <- fst %>%
  dplyr::mutate(module = factor(module, levels = mod_order)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, fst)
dxy_df <- dxy %>%
  dplyr::mutate(module = factor(module, levels = mod_order)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, dxy)

fst_boot <- fst_df %>% bootstrap_control_match(fst)
dxy_boot <- dxy_df %>% bootstrap_control_match(dxy)

# ---- new tiny core used by the wrappers below ----
if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a)) a else b

.mk_violin <- function(dfp, y, ylab, title, file = NULL,
                       star_q = 0.98, star_pad = 2,
                       y_range = NULL, fill_vals = NULL,
                       star_size = 7,
                       star_y = NULL   # <— fixed y for ALL stars (aligns LFMM & PCAdapt)
) {
  dfp <- dfp %>%
    dplyr::mutate(module_disp = module_to_disp(module)) %>%
    tidyr::drop_na({{ y }})
  
  mods_present <- intersect(mod_order_disp, unique(dfp$module_disp))
  dfp <- dfp %>% dplyr::mutate(module_disp = factor(module_disp, levels = mods_present))
  
  ctrl_med <- dfp %>%
    dplyr::filter(module_disp == "Control") %>%
    dplyr::summarise(med = stats::median({{ y }}, na.rm = TRUE)) %>% dplyr::pull(med)
  
  # p-values vs Control → stars
  df_for_test <- dfp %>% dplyr::transmute(module = module_disp, value = {{ y }})
  ptab <- wilcox_against_control(df_for_test, value) %>%
    dplyr::mutate(
      star = dplyr::case_when(
        is.na(p_raw) ~ "",
        p_raw < 0.001 ~ "***",
        p_raw < 0.01  ~ "**",
        p_raw < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::filter(star != "")
  
  # Fixed global y for star placement (if not provided, compute a safe one)
  if (is.null(star_y)) {
    star_y <- stats::quantile(dfp |> dplyr::pull({{ y }}), star_q, na.rm = TRUE) * star_pad
  }
  
  stars_df <- ptab %>%
    dplyr::transmute(
      module_disp = factor(group, levels = levels(dfp$module_disp)),
      label = star,
      ypos  = star_y
    )
  
  p <- ggplot(dfp, aes(x = module_disp, y = {{ y }}, fill = module_disp)) +
    geom_violin(trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    { if (!is.na(ctrl_med)) geom_hline(yintercept = ctrl_med, linetype = "dashed") } +
    { if (nrow(stars_df) > 0)
      geom_text(data = stars_df,
                aes(x = module_disp, y = ypos, label = label),
                inherit.aes = FALSE, fontface = "bold", size = star_size) } +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = (fill_vals %||% pal[mods_present]), guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1))
  
  # Ensure stars are visible
  if (!is.null(y_range)) {
    p <- p + coord_cartesian(ylim = y_range)
  } else if (nrow(stars_df)) {
    p <- p + expand_limits(y = star_y * 3)
  }
  
  if (!is.null(file)) {
    ggplot2::ggsave(filename = file, plot = p, width = 9, height = 5.5, dpi = 300)
  }
  invisible(list(plot = p, pvals = ptab))
}

# pick a fixed y for stars that’s above the highest box (adjust as needed)
fixed_star_y <- 0.23
plot_fst <- function(df, title = " ", ylab = expression(F[ST]),
                     file = NULL, y_range = NULL, star_q = 0.98, star_pad = 2) {
  .mk_violin(df, fst, ylab, title, file, star_q, star_pad, y_range,
             star_size = 7, star_y = fixed_star_y)
}

plot_dxy  <- function(df, title = " ", ylab = expression(d[xy]),
                      file = NULL, y_range = NULL, star_q = 0.98, star_pad = 1.05) {
  .mk_violin(df, dxy,  ylab, title, file, star_q, star_pad, y_range)
}

p_fst <- plot_fst(fst_boot, title = " ",
                  file = "fst_violin_box_boot.png",
                  y_range = c(0, 0.25))  # set e.g. c(0, 0.25) if you want
p_dxy <- plot_dxy(dxy_boot, title = " ",
                  file = "dxy_violin_box_boot.png",
                  y_range = c(0, 0.007))  # set e.g. c(0, 0.012)

fd_combined <- cowplot::plot_grid(p_fst$plot, p_dxy$plot, labels = c("A","B"), ncol = 2, align = "hv")
ggplot2::ggsave("fst_dxy_cowplot.png", fd_combined, width = 12, height = 5, dpi = 300)

cat("\n[Wilcoxon vs Control — FST]\n"); print(wilcox_against_control(fst_df, fst, reps = 1000))
cat("\n[Wilcoxon vs Control — DXY]\n"); print(wilcox_against_control(dxy_df, dxy, reps = 1000))

#-----------------------------------------------------------------------------------------
# ===================== ANGSD: θπ / θW / Fay’s H =====================
coast_FayH <- read_tsv("angsd_theta_by_gene_coast.tsv", show_col_types = FALSE)
north_FayH <- read_tsv("angsd_theta_by_gene_north.tsv", show_col_types = FALSE)

normalize_module <- function(x) tolower(sub("_genes$", "", as.character(x)))
mod_order <- c("control","egene","eqtl","lfmm","pcadapt","shared")

clean_angsd_by_module <- function(df){
  df %>%
    dplyr::rename(
      module = dplyr::any_of(c("module","Module")),
      gene   = dplyr::any_of(c("gene","Gene","id")),
      tW     = dplyr::any_of(c("tW","thetaW","ThetaW","watterson_theta","watterson")),
      tP     = dplyr::any_of(c("tP","thetaPi","ThetaPi","pi","theta_pi")),
      tH     = dplyr::any_of(c("tH","thetaH","ThetaH")),
      nSites = dplyr::any_of(c("nSites","nsites","Nsites","n_sites"))
    ) %>%
    dplyr::mutate(
      nSites = as.numeric(nSites),
      tW     = as.numeric(tW),
      tP     = as.numeric(tP),
      tH     = as.numeric(tH)
    ) %>%
    dplyr::filter(is.finite(nSites), nSites > 0) %>%
    dplyr::mutate(
      thetaW_per_site  = tW / nSites,
      thetaPi_per_site = tP / nSites,
      FayH_per_site    = (tP - tH) / nSites,
      module           = factor(normalize_module(module), levels = mod_order)
    )
}

coast_angsd <- clean_angsd_by_module(coast_FayH)
north_angsd <- clean_angsd_by_module(north_FayH)

coast_piA <- coast_angsd %>% dplyr::select(module, gene, thetaPi_per_site) %>% dplyr::rename(theta_pi = thetaPi_per_site)
coast_wA  <- coast_angsd %>% dplyr::select(module, gene, thetaW_per_site)  %>% dplyr::rename(theta_w  = thetaW_per_site)
coast_hA  <- coast_angsd %>% dplyr::select(module, gene, FayH_per_site)    %>% dplyr::rename(fayH     = FayH_per_site)
north_piA <- north_angsd %>% dplyr::select(module, gene, thetaPi_per_site) %>% dplyr::rename(theta_pi = thetaPi_per_site)
north_wA  <- north_angsd %>% dplyr::select(module, gene, thetaW_per_site)  %>% dplyr::rename(theta_w  = thetaW_per_site)
north_hA  <- north_angsd %>% dplyr::select(module, gene, FayH_per_site)    %>% dplyr::rename(fayH     = FayH_per_site)

coast_piA_boot  <- bootstrap_control_match(coast_piA,  theta_pi)
coast_wA_boot   <- bootstrap_control_match(coast_wA,   theta_w)
coast_hA_boot   <- bootstrap_control_match(coast_hA,   fayH)
north_piA_boot  <- bootstrap_control_match(north_piA,  theta_pi)
north_wA_boot   <- bootstrap_control_match(north_wA,   theta_w)
north_hA_boot   <- bootstrap_control_match(north_hA,   fayH)

ylab_theta_pi <- expression(theta[pi])
ylab_theta_w  <- expression(theta[W])
ylab_fayH     <- expression(Fay[H])

p_coast_pi_angsd <- plot_violin_box(coast_piA_boot,  theta_pi, "Coast", ylab_theta_pi, "angsd_coast_thetaPi_per_site.png")
p_coast_w_angsd  <- plot_violin_box(coast_wA_boot,   theta_w,  "Coast", ylab_theta_w,  "angsd_coast_thetaW_per_site.png")
p_north_pi_angsd <- plot_violin_box(north_piA_boot,  theta_pi, "North", ylab_theta_pi, "angsd_north_thetaPi_per_site.png")
p_north_w_angsd  <- plot_violin_box(north_wA_boot,   theta_w,  "North", ylab_theta_w,  "angsd_north_thetaW_per_site.png")

# ---- separate FayH plotter so you can coordinate-cartesian easily ----
plot_fayh <- function(df, title, ylab = expression(Fay[H]),
                      file = NULL, y_range = NULL, star_q = 0.995, star_pad = 1.20) {
  .mk_violin(df, fayH, ylab, title, file, star_q, star_pad, y_range)
}

p_coast_h_angsd <- plot_fayh(coast_hA_boot, "Coast",
                             file = "angsd_coast_fayH_per_site.png",
                             y_range = c(-0.02, 0.01))   # set e.g. c(-0.03, 0.01)
p_north_h_angsd <- plot_fayh(north_hA_boot, "North",
                             file = "angsd_north_fayH_per_site.png",
                             y_range = c(-0.02, 0.01))

angsd_grid <- cowplot::plot_grid(
  p_coast_pi_angsd$plot, p_north_pi_angsd$plot,
  p_coast_w_angsd$plot,  p_north_w_angsd$plot,
  labels = c("A","B","C","D"), label_size = 18, ncol = 2, align = "hv"
)
ggplot2::ggsave("angsd_theta_pi_coast_north_cowplot.png", angsd_grid, width = 12, height = 15, dpi = 300)

cat("\n[Wilcoxon — ANGSD θπ (Coast)]\n");  print(wilcox_against_control(coast_piA, theta_pi, reps = 1000))
cat("\n[Wilcoxon — ANGSD θW (Coast)]\n");  print(wilcox_against_control(coast_wA,  theta_w,  reps = 1000))
cat("\n[Wilcoxon — ANGSD FayH (Coast)]\n"); print(wilcox_against_control(coast_hA,  fayH,     reps = 1000))
cat("\n[Wilcoxon — ANGSD θπ (North)]\n");  print(wilcox_against_control(north_piA, theta_pi, reps = 1000))
cat("\n[Wilcoxon — ANGSD θW (North)]\n");  print(wilcox_against_control(north_wA,  theta_w,  reps = 1000))
cat("\n[Wilcoxon — ANGSD FayH (North)]\n"); print(wilcox_against_control(north_hA,  fayH,     reps = 1000))

#----------------------------------------------------------------------------------------
# ===================== CLR =====================
norm_module_key <- function(x){ tolower(gsub("_genes$", "", x)) }

plot_lr_by_module <- function(df, title, outfile = NULL,
                              # controls for star placement/size
                              star_size = 7,
                              star_y = NULL,          # if NULL, auto = 10% above max
                              star_q = 0.98,          # only used when star_y is NULL
                              star_y_pad = 1.10) {    # only used when star_y is NULL
  # prepare display labels and factor order
  dfp <- df %>%
    dplyr::mutate(module_disp = module_to_disp(module)) %>%
    tidyr::drop_na(LR)
  mods_present <- intersect(mod_order_disp, unique(dfp$module_disp))
  dfp <- dfp %>% dplyr::mutate(module_disp = factor(module_disp, levels = mods_present))
  
  # control median (dashed reference)
  ctl_med <- dfp %>%
    dplyr::filter(module_disp == "Control") %>%
    dplyr::summarise(med = median(LR, na.rm = TRUE)) %>%
    dplyr::pull(med)
  
  # Wilcoxon vs Control -> stars (use df *after* cleaning)
  df_for_test <- dfp %>% dplyr::transmute(module = module_disp, LR = LR)
  ptab <- wilcox_against_control(df_for_test, LR) %>%
    dplyr::mutate(
      star = dplyr::case_when(
        is.na(p_raw) ~ "",
        p_raw < 0.001 ~ "***",
        p_raw < 0.01  ~ "**",
        p_raw < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::filter(star != "")
  
  # one global y for all stars (aligns positions)
  if (is.null(star_y)) {
    # place stars slightly above the overall upper tail of LR
    ymax_q <- stats::quantile(dfp$LR, probs = star_q, na.rm = TRUE)
    star_y <- as.numeric(ymax_q) * star_y_pad
  }
  
  stars_df <- ptab %>%
    dplyr::transmute(
      module_disp = factor(group, levels = levels(dfp$module_disp)),
      label = star,
      y = star_y
    )
  
  p <- ggplot(dfp, aes(x = module_disp, y = LR, fill = module_disp)) +
    geom_violin(trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    { if (!is.na(ctl_med)) geom_hline(yintercept = ctl_med, linetype = "dashed") } +
    { if (nrow(stars_df) > 0)
      geom_text(data = stars_df, aes(x = module_disp, y = y, label = label),
                inherit.aes = FALSE, fontface = "bold", size = star_size) } +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = pal[mods_present], drop = FALSE, guide = "none") +
    labs(x = NULL, y = "CLR", title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    # ensure stars are visible
    expand_limits(y = star_y * 1.05)
  
  if (!is.null(outfile)) ggsave(outfile, p, width = 9, height = 5.5, dpi = 300)
  invisible(list(plot = p, pvals = ptab))
}


clr_data <- read_tsv(
  "/project/pi_brook_moyers_umb_edu/SF2/clr_analysis4/sf2_all_modules_by_gene_out/sf2_clr_max_by_gene.tsv",
  col_names = FALSE, show_col_types = FALSE
)
colnames(clr_data) <- c("module", "pop", "gene", "pos", "LR", "file")
clr_clean <- clr_data %>% dplyr::mutate(module_key = norm_module_key(module))
coast_lr <- clr_clean %>% dplyr::filter(pop == "coast") %>% dplyr::transmute(module = module_key, gene, LR)
north_lr <- clr_clean %>% dplyr::filter(pop == "north")  %>% dplyr::transmute(module = module_key, gene, LR)
coast_lr_boot <- bootstrap_control_match(coast_lr, LR)
north_lr_boot <- bootstrap_control_match(north_lr, LR)

p_coast_lr <- plot_lr_by_module(coast_lr_boot, "Coast", "coast_clr_violin_box.png")
p_north_lr <- plot_lr_by_module(north_lr_boot, "North", "north_clr_violin_box.png")

clr_grid <- cowplot::plot_grid(p_coast_lr$plot, p_north_lr$plot,
                               labels = c("A","B"), label_size = 18, ncol = 2, align = "hv")
ggplot2::ggsave("clr_violin_box_coast_north.png", clr_grid, width = 12, height = 5.5, dpi = 300)

# Big panel like the preview
clr_fay_grid <- cowplot::plot_grid(
  p_fst$plot, p_dxy$plot,
  p_coast_h_angsd$plot, p_north_h_angsd$plot,
  p_coast_lr$plot,  p_north_lr$plot,
  labels = c("A","B","C","D","E","F"), label_size = 18,
  ncol = 2, align = "hv"
)
ggplot2::ggsave("outliers_fst_dxy_fayh_clr_violin_box_coast_north.png", clr_fay_grid, width = 14, height = 12, dpi = 300)

cat("\n[Wilcoxon vs Control — CLR (Coast)]\n"); print(wilcox_against_control(coast_lr, LR, reps=1000))
cat("\n[Wilcoxon vs Control — CLR (North)]\n"); print(wilcox_against_control(north_lr,  LR, reps=1000))

# ===================== Connectivity =====================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})

# Safe loader that guarantees numeric 'connectivity'
read_connect <- function(path) {
  # readr is stricter & faster; parse_number strips commas/units if any
  readr::read_table2(
    file = path,
    col_names = c("gene","connectivity"),
    col_types = readr::cols(
      gene = readr::col_character(),
      connectivity = readr::col_character()
    )
  ) %>%
    mutate(
      connectivity = readr::parse_number(connectivity),
      gene = as.character(gene)
    )
}

# Wrapper around your shared .mk_violin core
plot_connectivity <- function(df,
                              title   = " ",
                              file    = NULL,
                              y_range = NULL,
                              star_q  = 0.98,
                              star_pad = 1.05) {
  .mk_violin(
    df,
    log_connectivity,
    ylab     = "Log10(Connectivity)",
    title    = title,
    file     = file,
    star_q   = star_q,
    star_pad = star_pad,
    y_range  = y_range
  )
}


lfmm_connectivity    <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_connectivity.txt")
pcadapt_connectivity <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/pcadapt_connectivity.txt")
eQTL_connectivity    <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/eQTL_connectivity.txt")
eGene_connectivity   <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/eGene_connectivity.txt")
control_connectivity <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/control_connectivity.txt")
shared_connectivity  <- read_connect("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/shared_connectivity.txt")

# If your files are clean numeric already, the simpler base read.table also works:
# read.table(..., col.names = c("gene","connectivity"), colClasses = c("character","numeric"))

# Build combined table
connect_df <- list(
  control = control_connectivity,
  lfmm    = lfmm_connectivity,
  pcadapt = pcadapt_connectivity,
  shared  = shared_connectivity,
  egene   = eGene_connectivity,
  eqtl    = eQTL_connectivity
) %>%
  dplyr::bind_rows(.id = "module") %>%
  { if (!exists("mod_order")) dplyr::mutate(., module = factor(module)) 
    else dplyr::mutate(., module = factor(module, levels = mod_order)) } %>%
  dplyr::mutate(
    connectivity = as.numeric(connectivity),      # ensure numeric
    log_connectivity = log10(connectivity + 1e-8) # add log-transformed column
  )

# Sanity checks (helpful if it fails again)
stopifnot(is.numeric(connect_df$log_connectivity))
# drop NA rows (wilcox.test can't handle NA)
connect_df <- dplyr::filter(connect_df, !is.na(log_connectivity))

# Optional: show any non-numeric leftovers before coercion (debug helper)
# bad <- readr::read_lines("/project/.../lfmm_connectivity.txt")
# dplyr::tibble(line = bad) %>% dplyr::filter(!str_detect(line, "^[^\\t]+\\t[-+]?[0-9]*\\.?[0-9]+([eE][-+]?\\d+)?$")) %>% print(n=50)

# Bootstrap Control to each module’s N
connect_boot <- connect_df %>% bootstrap_control_match(log_connectivity)

# Your wrapper should now work
p_connect <- plot_connectivity(
  connect_boot,
  title = " ",
  file = "connectivity_violin_box_boot.png",
  y_range = NULL
)

ggplot2::ggsave("connectivity.png", p_connect$plot, width = 10, height = 6, dpi = 300)


# ──────────────────────────────────────────────────────────────────────────────
# HaFT1 REGION (multiple genes): line+dot plots for π, θW, FST, DXY (windows)
# and Fay's H, CLR (per gene) — Coast vs North
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(data.table))

# Read region genes (chrom/start/end/gene) and prep positions
roi_genes <- read.csv("haft1_region_genes.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename(chromosome = chrom) %>%
  dplyr::mutate(
    chromosome = as.character(chromosome),
    start = as.numeric(start),
    end   = as.numeric(end)
  )
roi_gene_ids <- unique(roi_genes$gene)
roi_gene_pos <- roi_genes %>%
  dplyr::transmute(gene, chromosome, mid = 0.5 * (start + end))

# ---------------- helper: generic window-overlap via foverlaps ----------------
get_roi_windows_any <- function(df_windows, roi_tbl, value_col) {
  # df_windows must have: chromosome, window_pos_1, window_pos_2, gene, and `value_col`
  stopifnot(all(c("chromosome","window_pos_1","window_pos_2") %in% names(df_windows)))
  stopifnot(value_col %in% names(df_windows))
  
  x <- as.data.table(df_windows) %>%
    .[, .(chromosome = as.character(chromosome),
          start = as.numeric(window_pos_1),
          end   = as.numeric(window_pos_2),
          gene,
          value = as.numeric(get(value_col)))]
  setkey(x, chromosome, start, end)
  
  y <- as.data.table(roi_tbl) %>%
    .[, .(chromosome = as.character(chromosome),
          start = as.numeric(start),
          end   = as.numeric(end))]
  setkey(y, chromosome, start, end)
  
  ov <- foverlaps(x, y, nomatch = 0L)
  if (nrow(ov) == 0) return(tibble::tibble())
  out <- unique(ov[, .(chromosome, mid = 0.5*(start + end), gene, value)])
  tibble::as_tibble(out)
}

# ---------------- π & θW (window-based) ----------------
coast_div  <- get_roi_windows_any(coast_pi, roi_genes, "pi")  %>% dplyr::mutate(pop = "Coast", metric = "pi")
coast_watt <- get_roi_windows_any(coast_pi, roi_genes, "watterson_theta") %>% dplyr::mutate(pop = "Coast", metric = "watterson_theta")
north_div  <- get_roi_windows_any(north_pi, roi_genes, "pi")  %>% dplyr::mutate(pop = "North", metric = "pi")
north_watt <- get_roi_windows_any(north_pi, roi_genes, "watterson_theta") %>% dplyr::mutate(pop = "North", metric = "watterson_theta")

div_long <- dplyr::bind_rows(coast_div, coast_watt, north_div, north_watt) %>% dplyr::filter(is.finite(value))

p_div_line <- div_long %>%
  dplyr::arrange(metric, pop, mid) %>%
  ggplot(aes(x = mid, y = value, color = pop, group = pop)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y",
             labeller = labeller(metric = c(pi = expression(pi),
                                            watterson_theta = expression(theta[W])))) +
  scale_color_manual(values = c(Coast = "brown3", North = "steelblue")) +
  labs(x = "bp", y = NULL, color = "Population") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "right"
  )

# ---------------- FST & DXY (window-based) ----------------
coast_fst <- get_roi_windows_any(fst_df, roi_genes, "fst") %>% dplyr::mutate(pop = "Coast", metric = "fst")
north_fst <- get_roi_windows_any(fst_df, roi_genes, "fst") %>% dplyr::mutate(pop = "North", metric = "fst")  # fst_df already across pops
coast_dxy <- get_roi_windows_any(dxy_df, roi_genes, "dxy") %>% dplyr::mutate(pop = "Coast", metric = "dxy")
north_dxy <- get_roi_windows_any(dxy_df, roi_genes, "dxy") %>% dplyr::mutate(pop = "North", metric = "dxy")  # same note as above

fst_dxy_long <- dplyr::bind_rows(coast_fst, north_fst, coast_dxy, north_dxy) %>% dplyr::filter(is.finite(value))

p_fstdxy_line <- fst_dxy_long %>%
  dplyr::arrange(metric, pop, mid) %>%
  ggplot(aes(x = mid, y = value, color = pop, group = pop)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y",
             labeller = labeller(metric = c(fst = expression(F[ST]), dxy = expression(d[xy])))) +
  scale_color_manual(values = c(Coast = "brown3", North = "steelblue")) +
  labs(x = "bp", y = NULL, color = "Population") +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "right"
  )

# ---------------- Fay’s H & CLR (gene-based; join gene positions) ----------------
# ---------- Combine Fay's H and CLR into one long table ----------
fay_long <- dplyr::bind_rows(
  coast_angsd %>% dplyr::filter(gene %in% roi_gene_ids) %>%
    dplyr::transmute(pop = "Coast", gene, value = FayH_per_site),
  north_angsd %>% dplyr::filter(gene %in% roi_gene_ids) %>%
    dplyr::transmute(pop = "North", gene, value = FayH_per_site)
) %>%
  dplyr::left_join(roi_gene_pos, by = "gene") %>%
  dplyr::filter(is.finite(value), is.finite(mid)) %>%
  dplyr::mutate(metric = "FayH")

clr_long <- dplyr::bind_rows(
  coast_lr %>% dplyr::filter(gene %in% roi_gene_ids) %>%
    dplyr::transmute(pop = "Coast", gene, value = LR),
  north_lr %>% dplyr::filter(gene %in% roi_gene_ids) %>%
    dplyr::transmute(pop = "North", gene, value = LR)
) %>%
  dplyr::left_join(roi_gene_pos, by = "gene") %>%
  dplyr::filter(is.finite(value), is.finite(mid)) %>%
  dplyr::mutate(metric = "CLR")

scan_long <- dplyr::bind_rows(fay_long, clr_long)

# separate data just for the FayH=0 line
hline_df <- data.frame(metric = "FayH", yint = 0)

# ---------- One plot: Fay's H & CLR ----------
p_fay_clr <- scan_long %>%
  dplyr::arrange(metric, pop, mid) %>%
  ggplot(aes(x = mid, y = value, color = pop, group = pop)) +
  # zero line only on the FayH facet
  geom_hline(data = hline_df, aes(yintercept = yint),
             inherit.aes = FALSE, linewidth = 0.3, color = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ metric, scales = "free_y",
    labeller = as_labeller(c(FayH = "Fay's H", CLR = "CLR"))
  ) +
  scale_color_manual(values = c(Coast = "brown3", North = "steelblue")) +
  labs(x = "Position (bp)", y = NULL, color = "Population") +
  theme_bw() +
  theme(
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey90", color = NA),
    strip.text         = element_text(size = 13),
    axis.title.x       = element_text(size = 14),
    axis.title.y       = element_text(size = 14),
    axis.text.x        = element_text(size = 14),
    axis.text.y        = element_text(size = 14),
    legend.text        = element_text(size = 14),
    legend.title       = element_text(size = 14),
    legend.position    = "right"
  )

# optional: save
# ggsave("HaFT1_REGION_FayH_CLR_oneplot.png", p_fay_clr, width = 10, height = 4.8, dpi = 300)

# ---------------- Save individual panels (optional) ----------------
# ggsave("HaFT1_region_diversity_line.png", p_div_line, width=10, height=4.5, dpi=300)
# ggsave("HaFT1_region_FST_DXY_line.png", p_fstdxy_line, width=10, height=4.5, dpi=300)
# ggsave("HaFT1_region_FayH_line.png", p_fay_line, width=10, height=3.6, dpi=300)
# ggsave("HaFT1_region_CLR_line.png",  p_clr_line,  width=10, height=3.6, dpi=300)

# ---------------- Assemble 2×2 figure ----------------
haft1_panel_lines <- cowplot::plot_grid(
  p_div_line,          # A: π/θW
  p_fstdxy_line,       # B: FST/DXY
  p_fay_clr,        # D: CLR
  labels = c("A","B","C"), label_size = 16,
  ncol = 1, align = "hv"
)
ggsave("HaFT1_REGION_lines_pi_thetaW_FST_DXY_FayH_CLR_coast_vs_north.png",
       haft1_panel_lines, width = 16, height = 16, dpi = 300)

