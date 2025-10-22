#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(readr)
  library(stringr)
  library(vroom)
})

# =============================================================================
# CONFIG
# =============================================================================
indir  <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/pixy_out_genes/"
outdir <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] getwd():", getwd(), "\n")
cat("[INFO] outdir :", outdir, "\n")
cat("[INFO] writable outdir? ", file.access(outdir, 2) == 0, "\n")

# =============================================================================
# FILE DISCOVERY
# =============================================================================
pi_files   <- list.files(indir, "^pixy_genes_.*_pi\\.txt$",              full.names = TRUE)
dxy_files  <- list.files(indir, "^pixy_genes_.*_dxy\\.txt$",             full.names = TRUE)
fst_files  <- list.files(indir, "^pixy_genes_.*_fst\\.txt$",             full.names = TRUE)
wat_files  <- list.files(indir, "^pixy_genes_.*_watterson_theta\\.txt$", full.names = TRUE)
taj_files  <- list.files(indir, "^pixy_genes_.*_tajima_d\\.txt$",        full.names = TRUE)

# =============================================================================
# HELPERS
# =============================================================================
read_pixy <- function(paths, stat, main_col, keep_cols = NULL) {
  if (length(paths) == 0) return(tibble::tibble())
  df <- vroom::vroom(paths, delim = "\t", show_col_types = FALSE)
  df <- df %>%
    dplyr::rename(
      pop          = dplyr::any_of(c("pop","population")),
      chromosome   = dplyr::any_of(c("chromosome","chrom")),
      window_pos_1 = dplyr::any_of(c("window_pos_1","start","win_start")),
      window_pos_2 = dplyr::any_of(c("window_pos_2","end","win_end"))
    )
  keep_cols <- unique(c("pop","chromosome","window_pos_1","window_pos_2", main_col, keep_cols))
  df %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    dplyr::rename(!!stat := dplyr::all_of(main_col))
}

extract_module <- function(paths, stat){
  gsub(paste0("^pixy_genes_(.*)_", stat, "\\.txt$"), "\\1", basename(paths))
}

read_pixy_std <- function(path, stat, main_col, keep_cols = NULL){
  df <- vroom::vroom(path, delim = "\t", show_col_types = FALSE)
  df <- df %>%
    dplyr::rename(
      pop          = dplyr::any_of(c("pop","population")),
      pop1         = dplyr::any_of(c("pop1","population_1","pop_1")),
      pop2         = dplyr::any_of(c("pop2","population_2","pop_2")),
      chromosome   = dplyr::any_of(c("chromosome","chrom")),
      window_pos_1 = dplyr::any_of(c("window_pos_1","start","win_start")),
      window_pos_2 = dplyr::any_of(c("window_pos_2","end","win_end"))
    )
  keep_cols <- unique(c("pop","pop1","pop2","chromosome","window_pos_1","window_pos_2", main_col, keep_cols))
  df <- dplyr::select(df, dplyr::any_of(keep_cols))
  dplyr::rename(df, !!stat := dplyr::all_of(main_col))
}

read_bind_stat <- function(paths, stat, main_col, keep_cols = NULL){
  if (length(paths) == 0) return(tibble())
  modules <- extract_module(paths, stat)
  map2(paths, modules, ~{
    x <- read_pixy_std(.x, stat, main_col, keep_cols)
    x$module <- .y
    x
  }) %>% list_rbind()
}

fix_within <- function(x){
  if (!nrow(x)) return(x)
  if (!"pop" %in% names(x)) x$pop <- NA_character_
  if (is.list(x$pop) || is.data.frame(x$pop)) x$pop <- as.character(unlist(x$pop))
  x %>%
    mutate(
      pop          = as.character(pop),
      module       = as.character(module),
      chromosome   = as.character(chromosome),
      window_pos_1 = suppressWarnings(as.integer(window_pos_1)),
      window_pos_2 = suppressWarnings(as.integer(window_pos_2))
    )
}

fix_pairwise <- function(x){
  if (!nrow(x)) return(x)
  if (!"pop1" %in% names(x) && any(c("population_1","pop_1") %in% names(x)))
    x <- dplyr::rename(x, pop1 = dplyr::any_of(c("population_1","pop_1")))
  if (!"pop2" %in% names(x) && any(c("population_2","pop_2") %in% names(x)))
    x <- dplyr::rename(x, pop2 = dplyr::any_of(c("population_2","pop_2")))
  if (!"chromosome" %in% names(x) && "chrom" %in% names(x))
    x <- dplyr::rename(x, chromosome = chrom)
  
  x <- x %>%
    mutate(
      pop1         = tolower(trimws(as.character(pop1))),
      pop2         = tolower(trimws(as.character(pop2))),
      module       = as.character(module),
      chromosome   = as.character(chromosome),
      window_pos_1 = suppressWarnings(as.integer(window_pos_1)),
      window_pos_2 = suppressWarnings(as.integer(window_pos_2))
    ) %>%
    filter((pop1 == "coast" & pop2 == "north") | (pop1 == "north" & pop2 == "coast")) %>%
    mutate(pair = "coast__vs__north")
  
  x
}

save_df <- function(x, name, dir = outdir) {
  if (inherits(x, "data.frame") && nrow(x) > 0) {
    path <- file.path(dir, paste0(name, ".txt"))
    readr::write_tsv(x, path)
    cat("[OK] wrote:", path, " (", nrow(x), " rows)\n", sep = "")
  } else {
    cat("[SKIP] ", name, ": empty or not a data.frame\n", sep = "")
  }
}

# =============================================================================
# FIRST BLOCK: read simple pixy tables + write as-is
# =============================================================================
df_pi  <- read_pixy(pi_files,  "pi",  "avg_pi",
                    keep_cols = c("total_diffs","total_comps","total_missing"))

df_dxy <- read_pixy(dxy_files, "dxy", "avg_dxy",
                    keep_cols = c("total_diffs","total_comps","total_missing"))

df_fst <- read_pixy(fst_files, "fst", "avg_wc_fst",
                    keep_cols = c("a","b","c","n_sites")) %>%
  { .x <- .
  if (!"fst" %in% names(.x)) .x <- dplyr::rename(.x, fst = dplyr::any_of(c("wc_fst","hudson_fst","fst")))
  .x
  }

df_wat <- read_pixy(wat_files, "watterson_theta", "avg_watterson_theta",
                    keep_cols = c("raw_watterson_theta","no_sites","no_var_sites","weighted_no_sites")) %>%
  { .x <- .
  .x <- dplyr::rename(.x, raw_watterson_theta = dplyr::any_of(c("raw_theta","raw_watterson_theta")))
  .x <- dplyr::rename(.x, no_sites            = dplyr::any_of(c("num_sites","no_sites")))
  .x <- dplyr::rename(.x, no_var_sites        = dplyr::any_of(c("num_var_sites","no_var_sites")))
  .x <- dplyr::rename(.x, weighted_no_sites   = dplyr::any_of(c("num_weighted_sites","weighted_no_sites")))
  .x
  }

df_taj <- read_pixy(taj_files, "tajima_d", "tajima_d",
                    keep_cols = c("raw_pi","watterson_theta","d_stdev","num_sites","no_sites")) %>%
  dplyr::rename(no_sites = dplyr::any_of(c("num_sites","no_sites")))

save_df(df_pi,  "pixy_genes_pi")
save_df(df_dxy, "pixy_genes_dxy")
save_df(df_fst, "pixy_genes_fst")
save_df(df_wat, "pixy_genes_watterson_theta")
save_df(df_taj, "pixy_genes_tajima_d")

# =============================================================================
# SECOND BLOCK: module-aware reading, normalization, per-pop joins, BED matching
# =============================================================================
df_pi  <- read_bind_stat(pi_files,  "pi",  "avg_pi",
                         keep_cols = c("total_diffs","total_comps","total_missing")) %>% fix_within()
df_wat <- read_bind_stat(wat_files, "watterson_theta", "avg_watterson_theta",
                         keep_cols = c("raw_watterson_theta","no_sites","no_var_sites","weighted_no_sites")) %>%
  { .x <- .
  .x <- dplyr::rename(.x, raw_watterson_theta = dplyr::any_of(c("raw_theta","raw_watterson_theta")))
  .x <- dplyr::rename(.x, no_sites            = dplyr::any_of(c("num_sites","no_sites")))
  .x <- dplyr::rename(.x, no_var_sites        = dplyr::any_of(c("num_var_sites","no_var_sites")))
  .x <- dplyr::rename(.x, weighted_no_sites   = dplyr::any_of(c("num_weighted_sites","weighted_no_sites")))
  .x
  } %>% fix_within()

df_fst <- read_bind_stat(fst_files, "fst", "avg_wc_fst",
                         keep_cols = c("a","b","c","n_sites")) %>%
  { if (!"fst" %in% names(.)) dplyr::rename(., fst = dplyr::any_of(c("wc_fst","hudson_fst","fst"))) else . } %>%
  fix_pairwise()

df_dxy <- read_bind_stat(dxy_files, "dxy", "avg_dxy",
                         keep_cols = c("total_diffs","total_comps","total_missing")) %>% fix_pairwise()

# Build Coast/North per-pop tables (module included in join then dropped)
keys <- c("module","chromosome","window_pos_1","window_pos_2")

coast_pi  <- df_pi  %>% filter(pop == "Coast") %>% dplyr::select(all_of(keys), pi)
coast_wat <- df_wat %>% filter(pop == "Coast") %>% dplyr::select(all_of(keys), watterson_theta)
pair_fst  <- df_fst %>% dplyr::select(all_of(keys), fst)
pair_dxy  <- df_dxy %>% dplyr::select(all_of(keys), dxy)

coast_all <- list(coast_pi, coast_wat, pair_fst, pair_dxy) %>%
  keep(~ inherits(., "data.frame") && nrow(.) > 0) %>%
  reduce(full_join, by = keys) %>%
  arrange(chromosome, window_pos_1, window_pos_2)
coast_all$module <- NULL

north_pi  <- df_pi  %>% filter(pop == "North") %>% dplyr::select(all_of(keys), pi)
north_wat <- df_wat %>% filter(pop == "North") %>% dplyr::select(all_of(keys), watterson_theta)

north_all <- list(north_pi, north_wat, pair_fst, pair_dxy) %>%
  keep(~ inherits(., "data.frame") && nrow(.) > 0) %>%
  reduce(full_join, by = keys) %>%
  arrange(chromosome, window_pos_1, window_pos_2)
north_all$module <- NULL

# =============================================================================
# BEDs: read, tag modules, align to pixy coordinates (start1-1, end unchanged)
# =============================================================================
bed_dir   <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/"
bed_files <- list.files(bed_dir, "_genes.coords\\.bed$", full.names = TRUE)

beds <- map_dfr(bed_files, function(p) {
  mod <- str_replace(basename(p), "_genes.coords\\.bed$", "")
  read_tsv(p, col_names = c("chrom","start0","end1","gene"), col_types = "ciic") %>%
    mutate(module = mod, start1 = start0 + 1) %>%
    dplyr::select(module, gene, chrom, start1, end1)
})

beds2 <- beds %>%
  mutate(
    chrom      = as.character(chrom),
    start1     = as.integer(start1),
    end1       = as.integer(end1),
    pixy_start = start1 - 1L,  # key adjustment
    pixy_end   = end1
  )

# =============================================================================
# JOIN: subset π & θ by module via BED coordinates (per pop)
# =============================================================================
north_df <- north_all %>%
  dplyr::select(-any_of("module")) %>%
  mutate(
    chromosome   = as.character(chromosome),
    window_pos_1 = as.integer(window_pos_1),
    window_pos_2 = as.integer(window_pos_2)
  )

coast_df <- coast_all %>%
  dplyr::select(-any_of("module")) %>%
  mutate(
    chromosome   = as.character(chromosome),
    window_pos_1 = as.integer(window_pos_1),
    window_pos_2 = as.integer(window_pos_2)
  )

north_subset <- north_df %>%
  inner_join(
    beds2 %>% dplyr::select(module, gene, chrom, pixy_start, pixy_end),
    by = dplyr::join_by(
      chromosome   == chrom,
      window_pos_1 == pixy_start,
      window_pos_2 == pixy_end
    )
  ) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta) %>%
  arrange(module, chromosome, window_pos_1)

coast_subset <- coast_df %>%
  inner_join(
    beds2 %>% dplyr::select(module, gene, chrom, pixy_start, pixy_end),
    by = dplyr::join_by(
      chromosome   == chrom,
      window_pos_1 == pixy_start,
      window_pos_2 == pixy_end
    )
  ) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, pi, watterson_theta) %>%
  arrange(module, chromosome, window_pos_1)

readr::write_tsv(north_subset, file.path(outdir, "north_pi_watterson_by_module.tsv"))
readr::write_tsv(coast_subset, file.path(outdir, "coast_pi_watterson_by_module.tsv"))
cat("[OK] wrote north/coast pi+watterson subset files\n")

# =============================================================================
# JOIN: subset FST & DXY by module via BED coordinates (no pop split)
# =============================================================================
fst_by_module <- df_fst %>%
  dplyr::select(-any_of("module")) %>%
  mutate(
    chromosome   = as.character(chromosome),
    window_pos_1 = as.integer(window_pos_1),
    window_pos_2 = as.integer(window_pos_2)
  ) %>%
  inner_join(
    beds2,
    by = dplyr::join_by(
      chromosome   == chrom,
      window_pos_1 == pixy_start,
      window_pos_2 == pixy_end
    )
  ) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, fst) %>%
  arrange(module, chromosome, window_pos_1)

dxy_by_module <- df_dxy %>%
  dplyr::select(-any_of("module")) %>%
  mutate(
    chromosome   = as.character(chromosome),
    window_pos_1 = as.integer(window_pos_1),
    window_pos_2 = as.integer(window_pos_2)
  ) %>%
  inner_join(
    beds2,
    by = dplyr::join_by(
      chromosome   == chrom,
      window_pos_1 == pixy_start,
      window_pos_2 == pixy_end
    )
  ) %>%
  dplyr::select(module, gene, chromosome, window_pos_1, window_pos_2, dxy) %>%
  arrange(module, chromosome, window_pos_1)

readr::write_tsv(fst_by_module, file.path(outdir, "fst_by_module_from_bed.tsv"))
readr::write_tsv(dxy_by_module, file.path(outdir, "dxy_by_module_from_bed.tsv"))
cat("[OK] wrote fst/dxy by-module subset files\n")

# =============================================================================
# DONE
# =============================================================================
# Get ANGSD output and CLR 

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

BASE <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/angsd_per_gene_out"

# Modules and pops to scan (adjust if needed)
modules <- c("control_genes","eGene_genes","eQTL_genes","lfmm_genes","pcadapt_genes","shared_genes")
pops    <- c("coast","north")

# Build list of (module, pop, dir) that actually exist
dir_tbl <- expand_grid(module = modules, pop = pops) %>%
  mutate(dir = file.path(BASE, paste0(module, ".", pop))) %>%
  filter(dir.exists(dir))

if (nrow(dir_tbl) == 0L) {
  stop("No module/pop directories found under: ", BASE)
}

# Helper: read one pestPG file and normalize
read_pestpg <- function(f, module, pop) {
  # gene id from filename: gene.<GENE>.theta.thetas.idx.pestPG
  gene <- sub("^gene\\.(.*?)\\.theta\\.thetas\\.idx\\.pestPG$", "\\1", basename(f))
  df <- suppressMessages(read_tsv(f, show_col_types = FALSE, progress = FALSE))
  # first column often has a weird name like "(indexStart,...)" → rename to idx_info
  if (ncol(df) > 0 && grepl("^\\(", names(df)[1])) {
    names(df)[1] <- "idx_info"
  }
  # Ensure expected columns; keep a safe subset if present
  keep <- intersect(
    c("idx_info","Chr","WinCenter","tW","tP","tF","tH","tL",
      "Tajima","fuf","fud","fayh","zeng","nSites"),
    names(df)
  )
  df <- df[, keep, drop = FALSE]
  df$module <- module
  df$pop    <- pop
  df$gene   <- gene
  df
}

# Walk all dirs, bind
all_pest <- map2_dfr(dir_tbl$dir, seq_len(nrow(dir_tbl)), function(d, i) {
  module <- dir_tbl$module[i]
  pop    <- dir_tbl$pop[i]
  files  <- list.files(d, pattern = "\\.theta\\.thetas\\.idx\\.pestPG$", full.names = TRUE)
  if (length(files) == 0) return(tibble())
  bind_rows(lapply(files, read_pestpg, module = module, pop = pop))
})


# Reorder columns nicely
ord <- c("module","pop","gene","Chr","WinCenter","tW","tP","tF","tH","tL",
        "Tajima","fuf","fud","fayh","zeng","nSites","idx_info")
ord <- intersect(ord, names(all_pest))
 all_pest <- all_pest %>% dplyr::select(all_of(ord))

# Write combined and per-pop outputs
out_dir <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision/"
outfile_all   <- file.path(out_dir, "angsd_theta_by_gene_by_population.tsv")
outfile_coast <- file.path(out_dir, "angsd_theta_by_gene_coast.tsv")
outfile_north <- file.path(out_dir, "angsd_theta_by_gene_north.tsv")

write_tsv(all_pest,   outfile_all)
write_tsv(filter(all_pest, pop == "coast"), outfile_coast)
write_tsv(filter(all_pest, pop == "north"), outfile_north)

message("[OK] Wrote:\n  ", outfile_all,
        "\n  ", outfile_coast,
        "\n  ", outfile_north)

