suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readr); library(stringr)
})

# ---------- helper ----------
# Normalize: "chr01", "Ha412HOChr01", "Ha412HOChr1" -> "Ha412HOChr1"
normalize_chr <- function(x){
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("^Ha412HOChr0*([0-9]+)$", "Ha412HOChr\\1", x)
}

# ---------- read inputs ----------
lfmm_out    <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/lfmm_outliers_coords.csv")
pcadapt_out <- read.csv("pcadapt_outliers_coords.csv")
shared_out  <- read.table("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_pcadapt_common_pos.txt")
names(shared_out) <- c("CHR","POS")

# LFMM p-values (no header)
lfmm_pos <- read.table(
  "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_qval.txt",
  header = FALSE, stringsAsFactors = FALSE
) %>% dplyr::select(-V3, -V10, -V11)
names(lfmm_pos) <- c("chr","pos","REF","ALT","PC1","PC2","PC3","PC4")

# Collapse PCs -> one p per SNP (min across PCs)
lfmm_min <- lfmm_pos %>%
  mutate(chr = normalize_chr(chr), pos = as.integer(pos)) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "QVAL") %>%
  group_by(chr, pos) %>%
  summarise(QVAL = suppressWarnings(min(QVAL, na.rm = TRUE)), .groups = "drop") %>%
  mutate(QVAL = ifelse(is.infinite(QVAL), NA_real_, QVAL))

# pcadapt p-values (no header)
pcadapt_pos <- read.table(
  "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/pcadapt_qvalues.txt",
  header = FALSE, stringsAsFactors = FALSE
) %>% dplyr::select(-V3, -V7, -V8)
names(pcadapt_pos) <- c("chr","pos","REF","ALT","QVAL")
pcadapt_pos <- pcadapt_pos %>% mutate(chr = normalize_chr(chr), pos = as.integer(pos))

# ---------- merge to snpspos (must have chrom, pos) ----------
stopifnot(all(c("chrom","pos") %in% names(snpspos)))
base_pos <- snpspos %>% mutate(chrom = normalize_chr(chrom), pos = as.integer(pos))

# Per-method merged tables (keep separate to avoid cross-highlighting)
res_lfmm <- base_pos %>%
  left_join(lfmm_min %>% rename(chrom = chr, p = QVAL), by = c("chrom","pos")) %>%
  mutate(lp = ifelse(is.finite(p), -log10(p), NA_real_), method = "LFMM")

res_pcadapt <- base_pos %>%
  left_join(pcadapt_pos %>% rename(chrom = chr, p = QVAL), by = c("chrom","pos")) %>%
  mutate(lp = ifelse(is.finite(p), -log10(p), NA_real_), method = "pcadapt")

# ---------- chromosome offsets / centers ----------
all_chroms <- paste0("Ha412HOChr", 1:17)

chr_sizes <- base_pos %>%
  group_by(chrom) %>%
  summarise(chr_len = max(pos, na.rm = TRUE), .groups = "drop") %>%
  right_join(tibble(chrom = all_chroms), by = "chrom") %>%
  arrange(match(chrom, all_chroms)) %>%
  mutate(chr_len = ifelse(is.finite(chr_len), as.numeric(chr_len), 0))

chr_offsets <- chr_sizes %>% transmute(chrom, offset = dplyr::lag(cumsum(chr_len), default = 0))
axis_df     <- chr_sizes %>% left_join(chr_offsets, by = "chrom") %>%
  mutate(center = offset + chr_len/2)
axis_df_plot <- axis_df %>% filter(chr_len > 0)

# ---------- assign alternating chromosome colors ----------
chr_color_map <- axis_df_plot %>%
  mutate(
    chr_num   = as.integer(gsub("Ha412HOChr","", chrom)),
    chr_color = ifelse(chr_num %% 2L == 1L, "chr_grey", "chr_yellow")
  ) %>%
  dplyr::select(chrom, chr_color)

# attach chr_color to per-method tables (they already have bp_cum)
res_lfmm    <- res_lfmm    %>% left_join(chr_color_map, by = "chrom")
res_pcadapt <- res_pcadapt %>% left_join(chr_color_map, by = "chrom")

# ---------- outlier keys (as before, already normalized) ----------
lfmm_key <- lfmm_out %>%
  transmute(chrom = normalize_chr(CHR), pos = as.integer(POS)) %>% distinct() %>%
  mutate(is_lfmm = TRUE)

pcadapt_key <- pcadapt_out %>%
  transmute(chrom = normalize_chr(CHR), pos = as.integer(POS)) %>% distinct() %>%
  mutate(is_pcadapt = TRUE)

shared_key <- shared_out %>%
  transmute(chrom = normalize_chr(CHR), pos = as.integer(POS)) %>% distinct() %>%
  mutate(is_shared = TRUE)

# ---------- annotate per method: classify points ----------
# add cumulative bp to each method table
res_lfmm <- res_lfmm %>%
  left_join(chr_offsets, by = "chrom") %>%
  mutate(bp_cum = pos + offset)


res_pcadapt <- res_pcadapt %>%
  left_join(chr_offsets, by = "chrom") %>%
  mutate(bp_cum = pos + offset)

lfmm_plot <- res_lfmm %>%
  left_join(lfmm_key,  by = c("chrom","pos")) %>%
  left_join(shared_key, by = c("chrom","pos")) %>%
  mutate(cat = dplyr::case_when(
    is_shared ~ "shared",
    is_lfmm   ~ "lfmm",
    TRUE      ~ "other"
  ))

pcadapt_plot <- res_pcadapt %>%
  left_join(pcadapt_key, by = c("chrom","pos")) %>%
  left_join(shared_key,  by = c("chrom","pos")) %>%
  mutate(cat = dplyr::case_when(
    is_shared  ~ "shared",
    is_pcadapt ~ "pcadapt",
    TRUE       ~ "other"
  ))

plot_df <- bind_rows(lfmm_plot, pcadapt_plot) %>%
  mutate(method   = factor(method, levels = c("LFMM","pcadapt")),
         cat      = factor(cat,     levels = c("shared","lfmm","pcadapt","other")),
         chr_color= factor(chr_color, levels = c("chr_grey","chr_yellow"))) %>%
  filter(is.finite(lp))

# ---------- color palettes ----------
cols_outliers <- c(shared = "#4682B4", lfmm = "#A52A2A", pcadapt = "aquamarine4")
cols_chr      <- c(chr_grey = "grey40", chr_yellow = "#F4B400")   # odd=grey, even=yellow

# ---------- make chromosome 1 and 17 touch plot edges ----------
x_min <- min(axis_df_plot$offset, na.rm = TRUE)
x_max <- max(axis_df_plot$offset + axis_df_plot$chr_len, na.rm = TRUE)

# ---------- plot ----------
# ---------- ticks centered per chromosome actually present ----------
ticks_df <- axis_df %>%
  dplyr::filter(chr_len > 0, chrom %in% unique(plot_df$chrom)) %>%
  dplyr::arrange(match(chrom, paste0("Ha412HOChr", 1:17))) %>%
  dplyr::transmute(
    breaks = offset + chr_len/2,
    labels = gsub("Ha412HOChr", "", chrom)
  )

# Safety: x-limits spanning the plotted chromosomes
x_min <- min(axis_df$offset[axis_df$chr_len > 0], na.rm = TRUE)
x_max <- max(axis_df$offset[axis_df$chr_len > 0] + axis_df$chr_len[axis_df$chr_len > 0], na.rm = TRUE)

# ---------- plot ----------
p <- ggplot() +
  # background non-outliers alternate by chromosome
  geom_point(
    data = dplyr::filter(plot_df, cat == "other"),
    aes(x = bp_cum, y = lp, color = chr_color),
    size = 0.55, alpha = 0.65, na.rm = TRUE
  ) +
  # outliers on top: color by outlier type; shape by method (legend hidden)
  geom_point(
    data = dplyr::filter(plot_df, cat != "other"),
    aes(x = bp_cum, y = lp, color = cat, shape = method),
    size = 0.9, alpha = 0.95, na.rm = TRUE
  ) +
  # Color scale holds BOTH sets of mappings; legend shows only outlier types
  scale_color_manual(
    values = c(cols_outliers, cols_chr),
    breaks = names(cols_outliers),
    name   = "Outlier"
  ) +
  # Hide the method legend but keep shapes (circle vs triangle)
  scale_shape_manual(values = c(LFMM = 16, pcadapt = 17), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.9) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    expand = c(0, 0),
    breaks = ticks_df$breaks,             # centers
    labels = ticks_df$labels              # 1..17
  ) +
  labs(x = "Chromosome", y = expression(-log[10](q))) +
  theme_bw() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(size = 14, vjust = 1),
    axis.text.y        = element_text(size = 14),
    legend.title       = element_text(size = 12),
    legend.text        = element_text(size = 11)
  )

p
ggsave("manhattan_outliers_and_chralt.png", p, width = 12, height = 6, dpi = 300)

# ---------- Chr 6-only Manhattan (130–150 Mb) ----------
chr_target <- "Ha412HOChr6"
x_window   <- c(130e6, 150e6)
thr_p      <- 0.05
thr_lp     <- -log10(thr_p)

# color points that pass the p-value threshold even if not in the outlier CSVs
# color_threshold_hits <- TRUE

# Build  from method p-value tables
lfmm_chr6 <- lfmm_min %>%
  dplyr::rename(chrom = chr, p = QVAL) %>%
  dplyr::filter(chrom == chr_target, pos >= x_window[1], pos <= x_window[2]) %>%
  dplyr::mutate(lp = ifelse(is.finite(p), -log10(p), NA_real_), method = "LFMM") %>%
  dplyr::left_join(shared_key, by = c("chrom","pos")) %>%
  dplyr::left_join(lfmm_key,   by = c("chrom","pos")) %>%
  dplyr::mutate(
    is_list = is_shared | is_lfmm,
    is_thr  = is.finite(lp) & lp >= thr_lp,
    cat     = dplyr::case_when(
      is_shared ~ "shared",
      is_lfmm   ~ "lfmm",
      #color_threshold_hits & is_thr ~ "lfmm",  # treat LFMM threshold hits as lfmm color
      TRUE      ~ "other"
    )
  )

pcad_chr6 <- pcadapt_pos %>%
  dplyr::rename(chrom = chr) %>%
  dplyr::filter(chrom == chr_target, pos >= x_window[1], pos <= x_window[2]) %>%
  dplyr::mutate(
    lp = ifelse(is.finite(QVAL), -log10(QVAL), NA_real_),
    p  = QVAL,
    method = "pcadapt"
  ) %>%
  dplyr::left_join(shared_key,  by = c("chrom","pos")) %>%
  dplyr::left_join(pcadapt_key, by = c("chrom","pos")) %>%
  dplyr::mutate(
    is_list = is_shared | is_pcadapt,
    is_thr  = is.finite(lp) & lp >= thr_lp,
    cat     = dplyr::case_when(
      is_shared  ~ "shared",
      is_pcadapt ~ "pcadapt",
      #color_threshold_hits & is_thr ~ "pcadapt",  # treat pcadapt threshold hits as pcadapt color
      TRUE       ~ "other"
    )
  )

chr6_df <- dplyr::bind_rows(
  lfmm_chr6 %>% dplyr::select(pos, lp, p, method, cat),
  pcad_chr6 %>% dplyr::select(pos, lp, p, method, cat)
) %>%
  dplyr::filter(is.finite(lp)) %>%
  dplyr::mutate(
    cat    = factor(cat, levels = c("shared","lfmm","pcadapt","other")),
    method = factor(method, levels = c("LFMM","pcadapt"))
  )

cat("[Chr6] LFMM colored points: ",
    sum(chr6_df$method=="LFMM" & chr6_df$cat!="other"), "\n")
cat("[Chr6] pcadapt colored points: ",
    sum(chr6_df$method=="pcadapt" & chr6_df$cat!="other"), "\n")


# ---- plot ----
p_chr6_facets <- ggplot(chr6_df, aes(x = pos, y = lp)) +
  geom_point(
    data = dplyr::filter(chr6_df, cat == "other"),
    color = "grey60", size = 1.6, alpha = 0.6, shape = 16
  ) +
  geom_point(
    data = dplyr::filter(chr6_df, cat != "other"),
    aes(color = cat),
    size = 2.2, alpha = 0.95, shape = 16
  ) +
  scale_color_manual(values = cols_outliers, breaks = names(cols_outliers), name = "Outlier") +
  geom_hline(yintercept = thr_lp, linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(limits = x_window, expand = c(0,0), labels = scales::comma) +
  labs(x = "Chromosome 6 position (bp)", y = expression(-log[10](q))) +
  #facet_grid(rows = vars(method), scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )
# ggsave("chr6_130-150Mb_manhattan_facets.png", p_chr6_facets, width = 10, height = 7, dpi = 300)

ggsave("chr6_130-150Mb_manhattan.png", p_chr6_win, width = 10, height = 5, dpi = 300)



#Heat map
suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggnewscale)   # separate scale for the outlier strip
})

# ---------------- helpers ----------------
normalize_chr <- function(x){
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("^Ha412HOChr0*([0-9]+)$", "Ha412HOChr\\1", x)
}

# GT -> dosage (0,1,2), NA for missing
gt_to_dosage <- function(gt){
  gt <- gsub("\\|", "/", gt)
  out <- rep(NA_real_, length(gt))
  out[gt %in% c("0/0")]       <- 0
  out[gt %in% c("0/1","1/0")] <- 1
  out[gt %in% c("1/1")]       <- 2
  out
}

# make metadata IDs match VCF sample columns (e.g. "btm10-5" -> "btm10")
normalize_sample <- function(x){
  x <- tolower(trimws(as.character(x)))
  x <- sub("[._-][0-9]+$", "", x)      # strip trailing -5 / _03 / .2
  x <- sub("[._-]rep[0-9]+$", "", x)   # strip -rep2 if present
  x
}

# ---------------- inputs ----------------
vcf_file   <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/chr6_130-150.vcf.gz"
chr_target <- "Ha412HOChr6"
win_start  <- 130e6
win_end    <- 150e6

# sample metadata -> columns: sample (normalized IDs), pop ("north"/"coast")
sample_meta <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/arg_native_poplatlong.csv") %>%
  dplyr::filter(X %in% c("Coast", "North")) %>%
  dplyr::select(pop, X) %>%
  `names<-`(c("sample_raw","pop")) %>%
  mutate(
    pop    = tolower(pop),                # keep lower-case internally for ordering
    sample = normalize_sample(sample_raw)
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  dplyr::select(sample, pop)

# outlier tables
lfmm_out    <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/lfmm_outliers_coords.csv")
pcadapt_out <- read.csv("pcadapt_outliers_coords.csv")
shared_out  <- read.table("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_pcadapt_common_pos.txt",
                          col.names = c("CHR","POS"))

# ---------------- read VCF & slice region ----------------
vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE) %>%
  transmute(CHROM = normalize_chr(CHROM),
            POS   = as.integer(POS),
            ID    = ifelse(ID == ".", paste0(CHROM, ":", POS), ID))

keep <- which(fix$CHROM == chr_target & fix$POS >= win_start & fix$POS <= win_end)
stopifnot(length(keep) > 0)
vcf_win <- vcf[keep, ]

# ---------------- genotype matrix -> long (per sample) ----------------
gt_raw   <- vcfR::extract.gt(vcf_win, element = "GT", as.numeric = FALSE)  # variants x samples
dosage_w <- apply(gt_raw, 2, gt_to_dosage)
fix_win  <- fix[keep, ] %>% arrange(POS)
dosage_w <- dosage_w[order(fix[keep, "POS"]), , drop = FALSE]

vcf_samples <- colnames(gt_raw)
sample_key  <- tibble(sample_col = vcf_samples, sample = normalize_sample(vcf_samples))

geno_long <- as_tibble(dosage_w, .name_repair = "minimal") %>%
  mutate(POS = fix_win$POS) %>%
  relocate(POS) %>%
  setNames(c("POS", vcf_samples)) %>%
  pivot_longer(-POS, names_to = "sample_col", values_to = "dosage") %>%
  left_join(sample_key,  by = "sample_col") %>%
  left_join(sample_meta, by = "sample")

# quick diag
matched_n <- length(intersect(sample_key$sample, sample_meta$sample))
cat(sprintf("[INFO] metadata matched %d / %d VCF samples (%.1f%%)\n",
            matched_n, length(vcf_samples), 100*matched_n/length(vcf_samples)))

# ---------------- OUTLIERS-ONLY selection ----------------
lfmm_key <- lfmm_out    %>% transmute(CHROM = normalize_chr(CHR), POS = as.integer(POS), type = "lfmm")
pcad_key <- pcadapt_out %>% transmute(CHROM = normalize_chr(CHR), POS = as.integer(POS), type = "pcadapt")
shar_key <- shared_out  %>% transmute(CHROM = normalize_chr(CHR), POS = as.integer(POS), type = "shared")

# priority: shared > lfmm > pcadapt
out_annot <- bind_rows(lfmm_key, pcad_key, shar_key) %>%
  filter(CHROM == chr_target, POS >= win_start, POS <= win_end) %>%
  mutate(type = factor(type, levels = c("shared","lfmm","pcadapt"))) %>%
  arrange(type) %>%
  distinct(POS, .keep_all = TRUE) %>%
  dplyr::select(POS, type)

outlier_pos <- sort(unique(out_annot$POS))

# long table restricted to outlier columns (keep missing for now)
geno_out <- geno_long %>%
  filter(POS %in% outlier_pos, !is.na(pop)) %>%
  mutate(
    geno_cat = case_when(
      is.na(dosage) ~ NA_character_,
      dosage == 0   ~ "0/0",
      dosage == 1   ~ "0/1",
      dosage == 2   ~ "1/1"
    ),
    geno_cat = factor(geno_cat, levels = c("0/0","0/1","1/1", NA))
  )

# -------- NEW: drop samples with >=1 missing genotype at these outlier positions --------
miss_per_sample <- geno_out %>%
  group_by(sample) %>%
  summarise(n_missing = sum(is.na(geno_cat)), .groups = "drop")

keep_samples <- miss_per_sample %>%
  filter(n_missing == 0) %>%
  pull(sample)

geno_out <- geno_out %>%
  filter(sample %in% keep_samples)

# dense x index (no gaps among kept outliers)
pos_to_index <- tibble(POS = outlier_pos) %>%
  mutate(col = row_number())

geno_out <- geno_out %>%
  inner_join(pos_to_index, by = "POS")

# ---------------- order rows: North first, then Coast ----------------
# (capitalize for display, but keep ordering North -> Coast)
geno_out <- geno_out %>%
  mutate(
    pop     = factor(pop, levels = c("north","coast")),
    pop_lab = ifelse(pop == "north", "North", "Coast")
  )

row_order <- geno_out %>%
  distinct(sample, pop) %>%
  arrange(pop, sample) %>%
  pull(sample)

geno_out <- geno_out %>%
  mutate(
    sample = factor(sample, levels = rev(row_order)),  # top row = first
    row    = as.integer(sample)
  )

n_rows <- nlevels(geno_out$sample)

# right-side population labels (block midpoints) – labeled "North" / "Coast"
row_map <- geno_out %>%
  distinct(sample, pop, pop_lab) %>%
  mutate(y = as.integer(factor(sample, levels = levels(geno_out$sample))))

pop_mids <- row_map %>%
  group_by(pop, pop_lab) %>%
  summarise(mid = mean(range(y)), .groups = "drop") %>%
  arrange(factor(pop, levels = c("north","coast")))

# top strip: which outlier type each column is
strip_df <- out_annot %>%
  inner_join(pos_to_index, by = "POS") %>%
  transmute(snp_index = col, type)

# position for demarcation line between North and Coast
n_north <- sum(row_map$pop == "north")
sep_y   <- n_rows - n_north + 0.5   # horizontal line just below the last North sample
y_top   <- n_rows + 1               # top strip row

# ---------------- plot ----------------
geno_cols <- c(`0/0`="#88D36C", `0/1`="#1E9447", `1/1`="#0D5A2A")
out_cols  <- c(shared="#4682B4", lfmm="#A52A2A", pcadapt="aquamarine4")

# x ticks (≤6)
n_pos <- nrow(pos_to_index)
x_break_idx <- if (n_pos > 0) unique(round(seq(1, n_pos, length.out = min(6, n_pos)))) else integer(0)

p_heat <- ggplot(geno_out, aes(x = col, y = row)) +
  # genotype heatmap
  geom_raster(aes(fill = geno_cat)) +
  scale_fill_manual(values = geno_cols, name = "Genotype", na.translate = FALSE) +
  ggnewscale::new_scale_fill() +
  # outlier strip (thin band above the top row)
  { if (nrow(strip_df) > 0)
    geom_tile(
      data = strip_df,
      aes(x = snp_index, y = y_top, fill = type),
      height = 0.9, width = 1, inherit.aes = FALSE
    )
  } +
  scale_fill_manual(values = out_cols, name = "Outlier type") +
  # demarcation line between North and Coast
  geom_hline(yintercept = sep_y, color = "white", linewidth = 1.1) +
  # axes
  scale_x_continuous(
    breaks = pos_to_index$col[x_break_idx],
    labels = scales::comma(pos_to_index$POS[x_break_idx]),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, y_top + 0.5),
    breaks = seq_len(n_rows),
    labels = levels(geno_out$sample),   # left: sample names
    expand = c(0, 0),
    sec.axis = dup_axis(                # right: population blocks "North", "Coast"
      name   = " ",
      breaks = pop_mids$mid,
      labels = pop_mids$pop_lab
    )
  ) +
  labs(
    x = ("Chromosome 6 position (bp)"),
    y = " ") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", colour = NA),
    axis.text.x       = element_text(size = 14),
    axis.text.y       = element_text(size = 14),
    axis.title.x      = element_text(size = 14),
    axis.title.y.right= element_text(margin = margin(l = 6)),
    plot.title        = element_text(size = 18, face = "bold"),
    legend.position   = "right"
  ) +
  theme_bw()

print(p_heat)
ggsave("chr6_130-150Mb_genotype_heatmap_OUTLIERS_ONLY_no-missing-samples.png", p_heat, width = 10, height = 7, dpi = 300)

combined <- cowplot::plot_grid(p, p_chr6_facets, labels = c("A","B"), ncol = 1, align = "hv")
ggplot2::ggsave("HAFT1.png", combined, width = 10, height = 10, dpi = 300)


