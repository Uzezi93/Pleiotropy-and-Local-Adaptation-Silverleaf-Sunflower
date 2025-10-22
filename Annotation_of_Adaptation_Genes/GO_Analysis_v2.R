setwd("~/Documents/sunflower/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(rtracklayer)
})

# ------------------------------ helpers ---------------------------------

# read GFF and return a de-prefixed character vector of gene IDs
read_gff_genes <- function(path) {
  g <- readGFF(path)
  ids <- g$gene %>% as.character() %>% sub("^gene:", "", .)
  unique(ids[!is.na(ids) & ids != ""])
}

# run ORA for a gene set against a TERM2GENE map, with a defined background;
# saves table + barplot + dotplot
run_ora <- function(genes, name, TERM2GENE, background = NULL,
                    outdir = "ORA_results", width = 7, height = 4) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # keep only genes that appear in the mapping
  genes_use <- intersect(genes, TERM2GENE[[2]])
  if (length(genes_use) < 5) {
    warning(sprintf("[%s] too few genes after mapping (%d). Skipping.", name, length(genes_use)))
    return(invisible(NULL))
  }
  
  # background: intersect with TERM2GENE universal dataset
  univ <- if (is.null(background)) NULL else intersect(background, TERM2GENE[[2]])
  if (!is.null(univ) && length(univ) < 10) {
    warning(sprintf("[%s] background too small after intersection; ignoring background.", name))
    univ <- NULL
  }
  
  eg <- enricher(genes_use,
                 TERM2GENE = TERM2GENE,
                 universe  = univ,
                 pAdjustMethod = "BH")
  
  # if no enrichment, still save an empty table for bookkeeping
  tbl <- if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) as.data.frame(eg) else tibble()
  write.table(tbl, file = file.path(outdir, sprintf("%s_over_representation.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (nrow(tbl) > 0) {
    p_bar <- barplot(eg) + ggtitle(sprintf("%s – GO Over-representation", name))
    ggsave(file.path(outdir, sprintf("%s_barplot.tiff", name)), p_bar, width = width, height = height, dpi = 300)
    
    p_dot <- dotplot(eg) + ggtitle(sprintf("%s – GO Over-representation", name))
    ggsave(file.path(outdir, sprintf("%s_dotplot.tiff", name)), p_dot, width = width, height = height, dpi = 300)
  } else {
    message(sprintf("[%s] No enriched terms (after multiple-testing correction).", name))
  }
  
  invisible(eg)
}

# ------------------------------ inputs ----------------------------------

# gene sets (IDs will be cleaned inside reader)
outlier_genes <- read_gff_genes("~/Downloads/outlier_exons.gff")
shared_genes  <- read_gff_genes("~/Downloads/shared_exons.gff")
blue_genes    <- read_gff_genes("~/Downloads/blue_exons.gff")
pink_genes    <- read_gff_genes("~/Downloads/pink_exons.gff")
red_genes     <- read_gff_genes("~/Downloads/red_exons.gff")
yellow_genes  <- read_gff_genes("~/Downloads/yellow_exons.gff")   # <--- added
green_genes   <- read_gff_genes("~/Downloads/green_exons.gff")    # <--- already present
eGene_genes   <- read_gff_genes("eGenes_exons.gff")
expressed_genes <- read_gff_genes("~/Downloads/expressed_exons.gff")  # background

# GO mapping (filter to expressed background first)
# Your GAF: make sure 'Symbol' contains the same IDs as your GFF 'gene' field (no "gene:" prefix).
all_go <- read.delim("~/Downloads/GCF_002127325.2_HanXRQr2.0-SUNRISE_gene_ontology.gaf",
                     skip = 8, header = TRUE, check.names = FALSE)
all_go <- all_go %>%
  mutate(Symbol = sub("^gene:", "", as.character(Symbol))) %>%
  filter(Symbol %in% expressed_genes) %>%
  filter(!is.na(GO_ID) & GO_ID != "")

# get GO term names and join
go_terms <- clusterProfiler::go2term(unique(all_go$GO_ID)) %>%
  dplyr::rename(GO_ID = go_id, Term = Term)
all_go <- left_join(all_go, go_terms, by = "GO_ID")

# TERM2GENE: (term, gene) pairs
term2gene <- all_go %>%
  dplyr::select(Term, Symbol) %>%
  distinct() %>% na.omit()

# ------------------------------ run ORA for all sets ---------------------

bg <- expressed_genes  # background universe = all expressed genes

run_ora(outlier_genes, "All_Outliers", term2gene, background = bg)
run_ora(shared_genes,  "Shared_Outliers", term2gene, background = bg)

run_ora(blue_genes,    "Blue_Module",   term2gene, background = bg)
run_ora(pink_genes,    "Pink_Module",   term2gene, background = bg)
run_ora(red_genes,     "Red_Module",    term2gene, background = bg)
run_ora(yellow_genes,  "Yellow_Module", term2gene, background = bg)  # <--- added
run_ora(green_genes,   "Green_Module",  term2gene, background = bg)  # <--- added

run_ora(eGene_genes,   "eGenes",        term2gene, background = bg)

# ------------------------------ combine per-group TSVs into one file ------------------------------
outdir <- "ORA_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# pick up all per-group tables we just wrote
tsvs <- list.files(outdir, pattern = "_over_representation\\.tsv$", full.names = TRUE)

# read robustly (skip truly empty files), add Group from filename, and bind
dfs <- lapply(tsvs, function(f) {
  if (file.info(f)$size == 0L) return(NULL)  # empty file (no header/rows)
  df <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$Group <- sub("_over_representation\\.tsv$", "", basename(f))
  df
})
dfs <- Filter(Negate(is.null), dfs)

# always write a master file (empty or joined)
master_out <- file.path(outdir, "ALL_groups_over_representation.tsv")
if (length(dfs) == 0) {
  # write an empty, header-only table for consistency
  empty <- data.frame(
    Group=character(), ID=character(), Description=character(),
    GeneRatio=character(), BgRatio=character(),
    pvalue=double(), p.adjust=double(), qvalue=double(),
    geneID=character(), Count=integer(), stringsAsFactors = FALSE
  )
  write.table(empty, master_out, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  all_df <- dplyr::bind_rows(dfs) |>
    dplyr::relocate(Group, .before = 1)
  write.table(all_df, master_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Combined ORA written to: ", master_out)
