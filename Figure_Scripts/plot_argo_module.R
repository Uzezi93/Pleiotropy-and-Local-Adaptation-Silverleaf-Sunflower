#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggpubr)
  library(cowplot)
})

setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision/")

# ──────────────────────────────────────────────────────────────────────────────
# Globals: display order, palette, name mapping
# ──────────────────────────────────────────────────────────────────────────────
mod_order_disp <- c("Control","Blue","Pink","Yellow","Red","Green")

pal <- c(
  "Control" = "darkgrey",
  "Blue"    = "steelblue",
  "Pink"    = "deeppink",
  "Yellow"  = "goldenrod2",
  "Red"     = "firebrick3",
  "Green"   = "seagreen3"
)

module_to_disp <- function(x) {
  x0 <- tolower(as.character(x))
  x0 <- gsub("_genes$", "", x0)
  dplyr::recode(
    x0,
    "control_pool_all" = "Control",
    "control"          = "Control",
    "blue"  = "Blue",
    "pink"  = "Pink",
    "yellow"= "Yellow",
    "red"   = "Red",
    "green" = "Green",
    .default = stringr::str_to_title(x0)
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Helpers: bootstrap Control size & Wilcoxon tables (vs Control)
# ──────────────────────────────────────────────────────────────────────────────
bootstrap_control_match <- function(df, value_col, control_name="Control", seed=123){
  set.seed(seed)
  df <- df %>% mutate(module = module_to_disp(module))
  df_ctl <- df %>% filter(module == control_name) %>% tidyr::drop_na({{value_col}})
  df_oth <- df %>% filter(module != control_name) %>% tidyr::drop_na({{value_col}})
  if (nrow(df_ctl) == 0 || nrow(df_oth) == 0) return(df %>% tidyr::drop_na({{value_col}}))
  
  n_per_mod <- df_oth %>% count(module, name = "n_mod")
  ctl_pool  <- df_ctl %>% dplyr::select(all_of(setdiff(names(df_ctl), "module")))
  template_cols <- names(df)
  
  ctl_boot <- n_per_mod %>%
    group_by(module) %>%
    group_modify(function(.x, .y){
      n_m <- .x$n_mod[1]
      ctl_pool %>% slice_sample(n = n_m, replace = TRUE)
    }) %>% ungroup() %>%
    mutate(module = control_name) %>%
    dplyr::select(all_of(template_cols))
  
  bind_rows(df_oth, ctl_boot)
}

wilcox_against_control <- function(df, value_col, control_name="Control", reps=1000, seed=123){
  set.seed(seed)
  df <- df %>% mutate(module_disp = module_to_disp(module)) %>% tidyr::drop_na({{value_col}})
  present <- intersect(mod_order_disp, unique(df$module_disp))
  tests   <- setdiff(present, control_name)
  
  empty_out <- tibble(group=character(), n_test=integer(), p_raw=double(),
                      p_boot_median=double(), p_boot_min=double(), p_boot_max=double())
  if (length(tests) == 0 || !(control_name %in% present)) return(empty_out)
  
  ctl_vals <- df %>% filter(module_disp == control_name) %>% pull({{value_col}})
  if (!length(ctl_vals)) return(empty_out)
  
  purrr::map_dfr(tests, function(tg){
    x <- df %>% filter(module_disp == tg) %>% pull({{value_col}})
    if (!length(x)) {
      return(tibble(group=tg, n_test=0, p_raw=NA_real_,
                    p_boot_median=NA_real_, p_boot_min=NA_real_, p_boot_max=NA_real_))
    }
    p_raw <- tryCatch(wilcox.test(ctl_vals, x, exact = FALSE)$p.value, error=function(e) NA_real_)
    ps <- replicate(reps, {
      xb <- sample(ctl_vals, size = length(x), replace = TRUE)
      suppressWarnings(wilcox.test(xb, x, exact = FALSE)$p.value)
    })
    tibble(group = tg, n_test = length(x), p_raw = p_raw,
           p_boot_median = median(ps, na.rm = TRUE),
           p_boot_min    = min(ps, na.rm = TRUE),
           p_boot_max    = max(ps, na.rm = TRUE))
  })
}

# ──────────────────────────────────────────────────────────────────────────────
# π & Watterson θ — Coast / North
# ──────────────────────────────────────────────────────────────────────────────
coast <- read_tsv("coast_pi_watterson_by_color.tsv", show_col_types = FALSE)
north <- read_tsv("north_pi_watterson_by_color.tsv", show_col_types = FALSE)

coast_pi <- coast %>%
  mutate(module = factor(module_to_disp(module), levels = mod_order_disp)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta)

north_pi <- north %>%
  mutate(module = factor(module_to_disp(module), levels = mod_order_disp)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta)

# Bootstrap control sizes
coast_pi_boot  <- bootstrap_control_match(coast_pi,  pi) %>%
  mutate(module = factor(module, levels = mod_order_disp))
coast_wat_boot <- bootstrap_control_match(coast_pi,  watterson_theta) %>%
  mutate(module = factor(module, levels = mod_order_disp))
north_pi_boot  <- bootstrap_control_match(north_pi,  pi) %>%
  mutate(module = factor(module, levels = mod_order_disp))
north_wat_boot <- bootstrap_control_match(north_pi,  watterson_theta) %>%
  mutate(module = factor(module, levels = mod_order_disp))

# Coast: π
ctrl_med_plot <- coast_pi_boot %>% filter(module=="Control") %>%
  summarise(med = median(pi, na.rm = TRUE)) %>% pull(med)
ymax <- quantile(coast_pi_boot$pi, 0.995, na.rm = TRUE) * 1.05

p_coast_pi <- ggplot(coast_pi_boot, aes(module, pi, fill = module)) +
  geom_violin(adjust = 3, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = ctrl_med_plot, linetype = "dashed", color = "grey30") +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     ref.group = "Control", label.y = ymax * 11, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax)) +
  scale_fill_manual(values = pal[levels(coast_pi_boot$module)], drop = FALSE) +
  scale_y_continuous(trans = "sqrt") +
  labs(title = "Coast", x = NULL, y = "pi") +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1))

# Coast: θW
ctrl_med_plot_w <- coast_wat_boot %>% filter(module=="Control") %>%
  summarise(med = median(watterson_theta, na.rm = TRUE)) %>% pull(med)
ymax_w <- quantile(coast_wat_boot$watterson_theta, 0.995, na.rm = TRUE) * 1.05

p_coast_w <- ggplot(coast_wat_boot, aes(module, watterson_theta, fill = module)) +
  geom_violin(adjust = 3, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = ctrl_med_plot_w, linetype = "dashed", color = "grey30") +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     ref.group = "Control", label.y = ymax_w * 15, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_w)) +
  scale_fill_manual(values = pal[levels(coast_wat_boot$module)], drop = FALSE) +
  scale_y_continuous(trans = "sqrt") +
  labs(title = "Coast", x = NULL, y = "watterson theta") +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1))

# North: π
ctrl_med_plot_npi <- north_pi_boot %>% filter(module=="Control") %>%
  summarise(med = median(pi, na.rm = TRUE)) %>% pull(med)
ymax_npi <- quantile(north_pi_boot$pi, 0.995, na.rm = TRUE) * 1.05

p_north_pi <- ggplot(north_pi_boot, aes(module, pi, fill = module)) +
  geom_violin(adjust = 3, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = ctrl_med_plot_npi, linetype = "dashed", color = "grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_npi * 13, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_npi)) +
  scale_fill_manual(values = pal[levels(north_pi_boot$module)], drop = FALSE) +
  scale_y_continuous(trans = "sqrt") +
  labs(title = "North", x = NULL, y = "pi") +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1))

# North: θW
ctrl_med_plot_nw <- north_wat_boot %>% filter(module=="Control") %>%
  summarise(med = median(watterson_theta, na.rm = TRUE)) %>% pull(med)
ymax_nw <- quantile(north_wat_boot$watterson_theta, 0.995, na.rm = TRUE) * 1.05

p_north_w <- ggplot(north_wat_boot, aes(module, watterson_theta, fill = module)) +
  geom_violin(adjust = 3, scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = ctrl_med_plot_nw, linetype = "dashed", color = "grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_nw * 16, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_nw)) +
  scale_fill_manual(values = pal[levels(north_wat_boot$module)], drop = FALSE) +
  scale_y_continuous(trans = "sqrt") +
  labs(title = "North", x = NULL, y = "watterson theta") +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1))

# Combined π / θW panel
combined <- cowplot::plot_grid(
  p_coast_pi, p_north_pi,
  p_coast_w,  p_north_w,
  labels = c("A","B","C","D"), ncol = 2, align = "hv"
)
ggsave("module_pi_watterson_coast_north_cowplot.png", combined, width = 12, height = 10, dpi = 300)

# Wilcoxon tables (console)
cat("\n[Wilcoxon vs Control — COAST π]\n");  print(wilcox_against_control(coast_pi,  pi))
cat("\n[Wilcoxon vs Control — COAST θW]\n"); print(wilcox_against_control(coast_pi,  watterson_theta))
cat("\n[Wilcoxon vs Control — NORTH π]\n");  print(wilcox_against_control(north_pi,  pi))
cat("\n[Wilcoxon vs Control — NORTH θW]\n"); print(wilcox_against_control(north_pi,  watterson_theta))

# ──────────────────────────────────────────────────────────────────────────────
# FST / DXY (same look; no sqrt transform)
# ──────────────────────────────────────────────────────────────────────────────
fst <- read_tsv("fst_by_color_from_bed.tsv", show_col_types = FALSE)
dxy <- read_tsv("dxy_by_color_from_bed.tsv", show_col_types = FALSE)

fst_df <- fst %>% mutate(module = factor(module_to_disp(module), levels = mod_order_disp)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, fst)
dxy_df <- dxy %>% mutate(module = factor(module_to_disp(module), levels = mod_order_disp)) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, dxy)

fst_boot <- fst_df %>% bootstrap_control_match(fst) %>% mutate(module=factor(module, levels=mod_order_disp))
dxy_boot <- dxy_df %>% bootstrap_control_match(dxy) %>% mutate(module=factor(module, levels=mod_order_disp))

# FST
ctrl_med_fst <- fst_boot %>% filter(module=="Control") %>% summarise(med=median(fst,na.rm=TRUE)) %>% pull(med)
ymax_fst <- quantile(fst_boot$fst, 0.995, na.rm=TRUE) * 1.05
p_fst <- ggplot(fst_boot, aes(module, fst, fill = module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_fst, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_fst * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_fst)) +
  scale_fill_manual(values = pal[levels(fst_boot$module)], drop = FALSE) +
  labs(title = " ", x = NULL, y = expression(F[ST])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# DXY
ctrl_med_dxy <- dxy_boot %>% filter(module=="Control") %>% summarise(med=median(dxy,na.rm=TRUE)) %>% pull(med)
ymax_dxy <- quantile(dxy_boot$dxy, 0.995, na.rm=TRUE) * 1.05
p_dxy <- ggplot(dxy_boot, aes(module, dxy, fill = module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_dxy, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_dxy * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_dxy)) +
  scale_fill_manual(values = pal[levels(dxy_boot$module)], drop = FALSE) +
  labs(title = " ", x = NULL, y = expression(d[xy])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

fd_combined <- cowplot::plot_grid(p_fst, p_dxy, labels = c("A","B"), ncol = 2, align = "hv")
ggsave("module_fst_dxy_cowplot.png", fd_combined, width = 12, height = 5, dpi = 300)

cat("\n[Wilcoxon vs Control — FST]\n"); print(wilcox_against_control(fst_df, fst))
cat("\n[Wilcoxon vs Control — DXY]\n"); print(wilcox_against_control(dxy_df, dxy))

# ──────────────────────────────────────────────────────────────────────────────
# ANGSD per-site θπ / θW / Fay's H — Coast & North
# ──────────────────────────────────────────────────────────────────────────────

coast_FayH <- read_tsv("angsd_theta_by_module_coast.tsv", show_col_types = FALSE)
north_FayH <- read_tsv("angsd_theta_by_module_north.tsv", show_col_types = FALSE)

clean_angsd_by_module <- function(df){
  df %>%
    rename(
      module = any_of(c("module","Module")),
      gene   = any_of(c("gene","Gene","id")),
      tW     = any_of(c("tW","thetaW","ThetaW","watterson_theta","watterson")),
      tP     = any_of(c("tP","thetaPi","ThetaPi","pi","theta_pi")),
      tH     = any_of(c("tH","thetaH","ThetaH")),
      nSites = any_of(c("nSites","nsites","Nsites","n_sites"))
    ) %>%
    mutate(across(c(nSites,tW,tP,tH), as.numeric)) %>%
    filter(is.finite(nSites), nSites > 0) %>%
    mutate(
      module           = factor(module_to_disp(module), levels = mod_order_disp),
      thetaW_per_site  = tW / nSites,
      thetaPi_per_site = tP / nSites,
      FayH_per_site    = (tP - tH) / nSites
    )
}

coast_angsd <- clean_angsd_by_module(coast_FayH)
north_angsd <- clean_angsd_by_module(north_FayH)

coast_piA <- coast_angsd %>% dplyr::select(module, gene, thetaPi_per_site) %>% rename(theta_pi = thetaPi_per_site)
coast_wA  <- coast_angsd %>% dplyr::select(module, gene, thetaW_per_site)  %>% rename(theta_w  = thetaW_per_site)
coast_hA  <- coast_angsd %>% dplyr::select(module, gene, FayH_per_site)    %>% rename(fayH     = FayH_per_site)
north_piA <- north_angsd %>% dplyr::select(module, gene, thetaPi_per_site) %>% rename(theta_pi = thetaPi_per_site)
north_wA  <- north_angsd %>% dplyr::select(module, gene, thetaW_per_site)  %>% rename(theta_w  = thetaW_per_site)
north_hA  <- north_angsd %>% dplyr::select(module, gene, FayH_per_site)    %>% rename(fayH     = FayH_per_site)

coast_piA_boot  <- bootstrap_control_match(coast_piA,  theta_pi) %>% mutate(module=factor(module, levels=mod_order_disp))
coast_wA_boot   <- bootstrap_control_match(coast_wA,   theta_w)  %>% mutate(module=factor(module, levels=mod_order_disp))
coast_hA_boot   <- bootstrap_control_match(coast_hA,   fayH)     %>% mutate(module=factor(module, levels=mod_order_disp))
north_piA_boot  <- bootstrap_control_match(north_piA,  theta_pi) %>% mutate(module=factor(module, levels=mod_order_disp))
north_wA_boot   <- bootstrap_control_match(north_wA,   theta_w)  %>% mutate(module=factor(module, levels=mod_order_disp))
north_hA_boot   <- bootstrap_control_match(north_hA,   fayH)     %>% mutate(module=factor(module, levels=mod_order_disp))

# Coast θπ
ctrl_med_c_pi <- coast_piA_boot %>% filter(module=="Control") %>% summarise(med=median(theta_pi,na.rm=TRUE)) %>% pull(med)
ymax_c_pi <- quantile(coast_piA_boot$theta_pi, 0.995, na.rm=TRUE) * 1.05
p_coast_pi_angsd <- ggplot(coast_piA_boot, aes(module, theta_pi, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_c_pi, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_c_pi * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_c_pi)) +
  scale_fill_manual(values = pal[levels(coast_piA_boot$module)], drop = FALSE) +
  labs(title="Coast", x=NULL, y=expression(theta[pi])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# Coast θW
ctrl_med_c_w <- coast_wA_boot %>% filter(module=="Control") %>% summarise(med=median(theta_w,na.rm=TRUE)) %>% pull(med)
ymax_c_w <- quantile(coast_wA_boot$theta_w, 0.995, na.rm=TRUE) * 1.05
p_coast_w_angsd <- ggplot(coast_wA_boot, aes(module, theta_w, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_c_w, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_c_w * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_c_w)) +
  scale_fill_manual(values = pal[levels(coast_wA_boot$module)], drop = FALSE) +
  labs(title="Coast", x=NULL, y=expression(theta[W])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# Coast FayH
ctrl_med_c_h <- coast_hA_boot %>% filter(module=="Control") %>% summarise(med=median(fayH,na.rm=TRUE)) %>% pull(med)
ymax_c_h <- quantile(coast_hA_boot$fayH, 0.995, na.rm=TRUE) * 1.10
p_coast_h_angsd <- ggplot(coast_hA_boot, aes(module, fayH, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_c_h, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_c_h * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_c_h)) +
  scale_fill_manual(values = pal[levels(coast_hA_boot$module)], drop = FALSE) +
  labs(title="Coast", x=NULL, y=expression(Fay[H])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# North θπ
ctrl_med_n_pi <- north_piA_boot %>% filter(module=="Control") %>% summarise(med=median(theta_pi,na.rm=TRUE)) %>% pull(med)
ymax_n_pi <- quantile(north_piA_boot$theta_pi, 0.995, na.rm=TRUE) * 1.05
p_north_pi_angsd <- ggplot(north_piA_boot, aes(module, theta_pi, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_n_pi, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_n_pi * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_n_pi)) +
  scale_fill_manual(values = pal[levels(north_piA_boot$module)], drop = FALSE) +
  labs(title="North", x=NULL, y=expression(theta[pi])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# North θW
ctrl_med_n_w <- north_wA_boot %>% filter(module=="Control") %>% summarise(med=median(theta_w,na.rm=TRUE)) %>% pull(med)
ymax_n_w <- quantile(north_wA_boot$theta_w, 0.995, na.rm=TRUE) * 1.05
p_north_w_angsd <- ggplot(north_wA_boot, aes(module, theta_w, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_n_w, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_n_w * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_n_w)) +
  scale_fill_manual(values = pal[levels(north_wA_boot$module)], drop = FALSE) +
  labs(title="North", x=NULL, y=expression(theta[W])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# North FayH
ctrl_med_n_h <- north_hA_boot %>% filter(module=="Control") %>% summarise(med=median(fayH,na.rm=TRUE)) %>% pull(med)
ymax_n_h <- quantile(north_hA_boot$fayH, 0.995, na.rm=TRUE) * 1.10
p_north_h_angsd <- ggplot(north_hA_boot, aes(module, fayH, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_n_h, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_n_h * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_n_h)) +
  scale_fill_manual(values = pal[levels(north_hA_boot$module)], drop = FALSE) +
  labs(title="North", x=NULL, y=expression(Fay[H])) +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

angsd_grid <- cowplot::plot_grid(
  p_coast_pi_angsd, p_north_pi_angsd,
  p_coast_w_angsd,  p_north_w_angsd,
  labels = c("A","B","C","D"), label_size = 18, ncol = 2, align = "hv"
)
ggsave("module_angsd_theta_pi_coast_north_cowplot.png", angsd_grid, width = 12, height = 15, dpi = 300)

cat("\n[Wilcoxon vs Control — ANGSD θπ (Coast)]\n");  print(wilcox_against_control(coast_piA, theta_pi))
cat("\n[Wilcoxon vs Control — ANGSD θW (Coast)]\n");  print(wilcox_against_control(coast_wA,  theta_w))
cat("\n[Wilcoxon vs Control — ANGSD FayH (Coast)]\n"); print(wilcox_against_control(coast_hA,  fayH))
cat("\n[Wilcoxon vs Control — ANGSD θπ (North)]\n");  print(wilcox_against_control(north_piA, theta_pi))
cat("\n[Wilcoxon vs Control — ANGSD θW (North)]\n");  print(wilcox_against_control(north_wA,  theta_w))
cat("\n[Wilcoxon vs Control — ANGSD FayH (North)]\n"); print(wilcox_against_control(north_hA,  fayH))

# ──────────────────────────────────────────────────────────────────────────────
# CLR — Coast & North (same look; no sqrt transform)
# ──────────────────────────────────────────────────────────────────────────────
norm_module_key <- function(x){ tolower(gsub("_genes$", "", x)) }

clr_data <- read_tsv(
  "/project/pi_brook_moyers_umb_edu/SF2/clr_analysis5/sf2_color_modules_by_gene_out/sf2_clr_max_by_gene.tsv",
  col_names = FALSE, show_col_types = FALSE
)
colnames(clr_data) <- c("module", "pop", "gene", "pos", "LR", "file")
clr_clean <- clr_data %>% mutate(module_key = norm_module_key(module))

coast_lr <- clr_clean %>% filter(pop == "coast") %>% transmute(module = module_key, gene, LR)
north_lr <- clr_clean %>% filter(pop == "north") %>% transmute(module = module_key, gene, LR)

coast_lr_boot <- bootstrap_control_match(coast_lr, LR) %>% mutate(module=factor(module_to_disp(module), levels=mod_order_disp))
north_lr_boot <- bootstrap_control_match(north_lr, LR) %>% mutate(module=factor(module_to_disp(module), levels=mod_order_disp))

# Coast CLR
ctrl_med_clr_c <- coast_lr_boot %>% filter(module=="Control") %>% summarise(med=median(LR,na.rm=TRUE)) %>% pull(med)
ymax_clr_c <- quantile(coast_lr_boot$LR, 0.995, na.rm=TRUE) * 1.05
p_coast_lr <- ggplot(coast_lr_boot, aes(module, LR, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_clr_c, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_clr_c * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_clr_c)) +
  scale_fill_manual(values = pal[levels(coast_lr_boot$module)], drop = FALSE) +
  labs(title="Coast", x=NULL, y="CLR") +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

# North CLR
ctrl_med_clr_n <- north_lr_boot %>% filter(module=="Control") %>% summarise(med=median(LR,na.rm=TRUE)) %>% pull(med)
ymax_clr_n <- quantile(north_lr_boot$LR, 0.995, na.rm=TRUE) * 1.05
p_north_lr <- ggplot(north_lr_boot, aes(module, LR, fill=module)) +
  geom_violin(adjust=3, scale="width", trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_hline(yintercept=ctrl_med_clr_n, linetype="dashed", color="grey30") +
  stat_compare_means(method="wilcox.test", label="p.signif", ref.group="Control",
                     label.y = ymax_clr_n * 0.92, hide.ns = TRUE, size =7) +
  coord_cartesian(ylim = c(NA, ymax_clr_n)) +
  scale_fill_manual(values = pal[levels(north_lr_boot$module)], drop = FALSE) +
  labs(title="North", x=NULL, y="CLR") +
  guides(fill="none") + theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35,hjust=1))

clr_grid <- cowplot::plot_grid(p_coast_lr, p_north_lr, labels = c("A","B"), ncol = 2, align = "hv")
ggsave("clr_violin_box_coast_north.png", clr_grid, width = 12, height = 5.5, dpi = 300)

cat("\n[Wilcoxon vs Control — CLR (Coast)]\n"); print(wilcox_against_control(coast_lr, LR))
cat("\n[Wilcoxon vs Control — CLR (North)]\n"); print(wilcox_against_control(north_lr,  LR))


# Big panel like the preview
clr_fay_grid <- cowplot::plot_grid(
  p_fst, p_dxy,
  p_coast_h_angsd, p_north_h_angsd,
  p_coast_lr,  p_north_lr,
  labels = c("A","B","C","D","E","F"), label_size = 18,
  ncol = 2, align = "hv"
)
ggsave("fst_dxy_fayh_clr_violin_box_coast_north.png", clr_fay_grid, width = 12, height = 10, dpi = 300)


