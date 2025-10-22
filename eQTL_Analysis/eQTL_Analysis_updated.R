---
  title: "eQTL_analysis"
author: "Uzezi Okinedo"
date: "09/4/2023"
output: html_document
---


setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/eQTL_analysis/")

# load tidyverse package
library(tidyverse)
library(tidyr)
library(MatrixEQTL)
library(trio)
library(broom)
library(Biobase)
library(qvalue)
library(tidyverse)
library(data.table)

# Process SNPs for eQTL analysis
snps1 <- read.table("~/Downloads/genos.012", sep = "\t", header = TRUE, row.names = NULL)
snps <- data.frame(t(snps1))

# Edit snp dataframe to reflect same number of rows as loci dataset
snps <- snps[-c(1, 113626, 113627), ]
snp_colname <- read.table("~/Downloads/genos.012.indv", sep = "\t", header = FALSE, skip = 1)
snp_loc <- read.table("~/Downloads/genos.012.pos", sep = "\t", header = FALSE, skip = 1)
snp_loc <- snp_loc[-c(1), ]
row.names(snp_loc) <- make.names(snp_loc$V1, unique = TRUE)
colnames(snps) <- snp_colname$V1
row.names(snps) <- row.names(snp_loc)
snps <- snps[, order(names(snps))]

# snp column to snp_loc data
snp_loc <- rownames_to_column(snp_loc, var = "RowNames")
colnames(snp_loc) <- c("snps", "chrom", "pos")

# Remove genebank accessions from snps dataset
snps$Ames449 <- NULL
snps$Ames695 <- NULL
snps$ARG1805<- NULL
snps$ARG1820 <- NULL
snps$ARG1834 <- NULL

# Remove missing genotypes
snps <- subset(snps, !rowSums(snps < 0)) # this seems to be biased against the reference
#k <- setdiff(snp_loc$snps, row.names(snps))
#snp_loc <- snp_loc[!snp_loc$snps %in% k, ]

# Repair snp colnames
colnames(snps) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(snps))
colnames(snps) <- gsub("-", ".", colnames(snps))
snps$btm7B.14 <- NULL
snps$arg11B.11 <- NULL
snps$btm5.1 <- NULL

# Compare snp data set to expression data set
expr <- read.table("~/Downloads/expression_vst_matrix.txt", sep = "\t", header = TRUE)

expr2 <- expr

#res_data <- read.table("expr_residual_matrix.txt")
# expr2 <- data.frame(t(expr))

# expr2 <- janitor::row_to_names(expr2, row_number = 1)

expr2$btm5.1 <- NULL
expr2$arg11B.11 <- NULL
expr2$btm9.4 <- NULL
expr2$btm7B.14 <- NULL
expr2$Ames449 <- NULL
expr2$Ames695 <- NULL
expr2$ARG1805<- NULL
expr2$ARG1820 <- NULL
expr2$ARG1834 <- NULL
#expr$X <- NULL

#row.names(expr2) <- expr$Gene

# Convert expr to numeric
# expr <- as.data.frame(lapply(expr, function(x) as.numeric(as.character(x))))

# Convert the row names of the original data frame into column names
# expr$Gene_id <- rownames(expr)

expr_cols <-names(expr2[, 2:19])

expr_cols <- gsub("-", ".", expr_cols)

snps <- snps %>% dplyr::select(all_of(expr_cols))


# Ensure eqtl is run using genes with in snp dataset
# get gene names in snp dataset
gene_names <- read.table("~/Downloads/variable_genes2.txt")


# 1. Remove the "gene:" prefix if present
expr2$Gene <- sub("^gene:", "", expr2$Gene)
rownames(expr2) <- expr2$Gene

# 2. Remove everything after the first dot, if it exists (e.g., version suffixes)
# expr$Gene_id_clean <- sub("\\..*$", "", expr$Gene_id_clean)

# Preview cleaned IDs
# head(expr$Gene_id_clean, 20)

## Subset expr to only rows whose Gene_id_clean are in gene_names$clean
expr_subset <- expr2[expr2$Gene %in% gene_names$V1, ]

## Check
nrow(expr_subset)
head(expr_subset$Gene)

#rownames(expr_subset) <- expr_subset$Gene_id_clean
expr_subset$Gene <- NULL
expr_subset$Gene <- NULL

# colnames(expr_subset)[1:18] <- gsub("-", ".", colnames(expr_subset)[1:18])

#expr_subset <- expr_subset %>% select(expr_cols)


# Remove missing locus?
gene_loc <- read.table("~/Downloads/gene_loc.txt", header = TRUE)

# --- Clean IDs in gene_loc ---
gene_loc$Geneid_clean <- gene_loc$Geneid

# Replace "gene." with "gene:" for consistency
gene_loc$Geneid_clean <- sub("^gene\\.", "gene:", gene_loc$Geneid_clean)
gene_loc$Geneid_clean <- sub("^gene:", "", gene_loc$Geneid_clean)

# --- Now both matchable ---
head(gene_loc$Geneid_clean)
head(rownames(expr_subset))

#rownames(gene_loc) <- gene_loc$Geneid
h <- setdiff(gene_loc$Geneid_clean, rownames(expr_subset))
gene_loc <- gene_loc[!gene_loc$Geneid_clean %in% h, ]
gene_loc$Geneid_clean <- NULL
gene_loc$Geneid <- sub("^gene\\.", "gene:", gene_loc$Geneid)
gene_loc$Geneid <- sub("^gene:", "", gene_loc$Geneid)

# Prepare covariate file
cov <- read.delim("~/Downloads/covariates_pc12.txt", sep = "\t", fill = TRUE)

# Make cov column names consistent with expr_cols
# names(cov) <- gsub("-", ".", names(cov))

# Now subset using expr_cols
cov <- cov %>%
  dplyr::select(all_of(expr_cols))


# Write snps, loci, and datasets into file
write.table(snps, file = "~/Documents/eQTL_snps.txt", sep = "\t", quote = FALSE)
write.table(snp_loc, file = "~/Documents/snp_loc.txt", sep = "\t", quote = FALSE)
write.table(expr_subset, file = "~/Documents/expr_residuals.txt", sep = "\t", quote = FALSE)
write.table(gene_loc, file = "~/Documents/gene_loc.txt", quote = FALSE)
write.table(cov, file = "~/Documents/covariates.txt", sep = "\t", quote = FALSE)

# Location of the package with the data files.
base.dir = find.package('MatrixEQTL')
base.dir = '~/Documents/'

# Settings
useModel = MatrixEQTL::modelLINEAR

# Genotype file name
SNP_file_name = paste(base.dir, "eQTL_snps.txt", sep = "")
snps_location_file_name = paste(base.dir, "snp_loc.txt", sep = "")

# Gene expression file name
expression_file_name = paste(base.dir, "expr_residuals.txt", sep = "")
gene_location_file_name = paste(base.dir, "gene_loc.txt", sep = "")

# Covariates file name
covariates_file_name = paste(base.dir, "covariates.txt", sep = "")

# Output file name
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-3
pvOutputThreshold_tra = 1e-3

# Error covariance matrix
errorCovariance = numeric()

# Distance for local gene-SNP pairs
cisDist = 1e6

# Load genotype data
snps = MatrixEQTL::SlicedData$new()
snps$fileDelimiter = "\t"       # the TAB character
snps$fileOmitCharacters = "NA"  # denote missing values;
snps$fileSkipRows = 1           # one row of column labels
snps$fileSkipColumns = 1        # one column of row labels
snps$fileSliceSize = 2000       # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

# Load gene expression data
gene = MatrixEQTL::SlicedData$new()
gene$fileDelimiter = "\t"       # the TAB character
gene$fileOmitCharacters = "NA"  # denote missing values;
gene$fileSkipRows = 1           # one row of column labels
gene$fileSkipColumns = 1        # one column of row labels
gene$fileSliceSize = 2000       # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

# Load covariates
cvrt = MatrixEQTL::SlicedData$new()
cvrt$fileDelimiter = "\t"       # the TAB character
cvrt$fileOmitCharacters = "NA"  # denote missing values;
cvrt$fileSkipRows = 1           # one row of column labels
cvrt$fileSkipColumns = 1        # one column of row labels
if(length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

# Remove SNPs with low MAF
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]]
  maf.list[[sl]] = rowMeans(slice, na.rm = TRUE) / 2
  maf.list[[sl]] = pmin(maf.list[[sl]], 1 - maf.list[[sl]])
}
maf = unlist(maf.list)

# Verify length of maf vector
cat('Length of MAF vector:', length(maf), '\n')
cat('Number of rows in SNPs:', nrow(snps), '\n')

# Check distribution of MAF
hist(maf)

# Filter maf values greater than 0.05
# filtered_maf <- maf[maf < 0.1]

# Plot histogram with appropriate breaks
# hist(filtered_maf, breaks = seq(0.1, max(filtered_maf), length.out = 100),
#     main = "Filtered MAF Histogram", xlab = "maf", col = "gray")

# hist(maf[maf < 0.05], seq(0, 0.05, length.out = 100))

# Filter SNPs by MAF
cat('SNPs before filtering:', nrow(snps), '\n')
# snps$RowReorder(maf > 0.05)
snps$RowReorderSimple(maf > 0.1)
cat('SNPs after filtering:', nrow(snps), '\n')

# Remove outliers
#for( sl in 1:length(gene) ) {
# mat = gene[[sl]];
# mat = t(apply(mat, 1, rank, ties.method = "average"));
# mat = qnorm(mat / (ncol(gene)+1));
# gene[[sl]] = mat;
# }
# rm(sl, mat);

# Load SNP and gene positions
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

# Run the analysis
me = MatrixEQTL::Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  #cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

unlink(output_file_name_tra)
unlink(output_file_name_cis)

# Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected local eQTLs:', '\n')
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n')
show(me$trans$eqtls)

trans_eQTL <- me$trans$eqtls
trans_eQTL2 <- trans_eQTL %>% filter(FDR < 0.05)

write.table(trans_eQTL2, file = "~/Downloads/transeqtls.txt", sep = "\t", quote = F, row.names = F)
# Get eQTLs and eGenes
# trans_eQTL <- read.table("/project/umb_brook_moyers/argo_network/trans_eQTL.txt",sep="\t",
#                                      header=T,row.names=NULL, fill=TRUE) 

# Get eGenes and eQTLs
eGenes <- data.frame(unique(trans_eQTL2$gene))
write.table(eGenes, file = "eGenes.txt", sep = "\t", quote = F, row.names = F)

eQTLs <- data.frame(unique(trans_eQTL2$snps))
write.table(eQTLs, file = "eQTLs.txt")

# Extract eQTL and eGenes positions
eQTL_positions <- snpspos %>%
  filter(snps %in% trans_eQTL2$snps) %>%
  distinct()

eQTL_positions$snps <- NULL

write.table(eQTL_positions, file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/eQTL_position2.txt", sep = "\t", quote = F, row.names = F)

eGenes_positions <- genepos %>%
  filter(Geneid %in% trans_eQTL2$gene) %>%
  distinct() 

eGenes_positions$Geneid <- NULL

write.table(eGenes_positions, file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/eGenes_position2.txt", sep = "\t", quote = F, row.names = F)

eGenes <- data.frame(trans_eQTL2$gene) 

eGenes <- eGenes %>%
  separate(trans_eQTL2.gene, into = c("Gene_ID", "ID", NA)) %>%
  distinct()

eGenes$Gene_ID <- NULL

# Genes with sufficient expression
suff_gene <- data.frame(row.names(expr_subset))

suff_gene <- suff_gene %>%
  #separate(row.names.expr_subset., into = c("NA", "Gene_ID")) %>%
  distinct()

write.table(suff_gene, file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/all_genes_eQTL.txt", sep = "\t", quote = FALSE)

length(unique(suff_gene$row.names.expr_subset.)) # 14757


# Genes with significant eQTL associations = 24
eGenes$Gene_ID <- NULL
write.table(eGenes, file = "eGenes.txt", sep = "\t", quote = F, row.names = F)


## Check if outliers on average have more eGenes than non-outliers
## ── Libraries ────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)

stopifnot(all(c("snps","chrom","pos") %in% names(snpspos)))

## User params
FDR_cut   <- 0.05
gff_path  <- "~/Downloads/HAN412_Eugene_curated_v1_1.gff3"   # <-- set your GFF path
lfmm_path    <- "~/Downloads/lfmm_pos.txt"
pcadapt_path <- "~/Downloads/pcadapt_pos.txt"

## - Functions ──────────────────────────────────────────────────────────────────
read_outlier_pairs <- function(path, zero_based = FALSE) {
  x <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(x) == 1) {
    x <- tidyr::separate(x, V1, into = c("chr","pos"), sep = "[:\\s]+", remove = TRUE, convert = TRUE)
  } else {
    names(x)[1:2] <- c("chr","pos")
  }
  x <- dplyr::mutate(x,
                     chr = as.character(chr),
                     pos = as.integer(pos) + ifelse(zero_based, 1L, 0L))
  dplyr::distinct(x, chr, pos)
}

normalize_chr <- function(x) {
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("^Ha412HOChr0?([0-9]+)$", "Ha412HOChr\\1", x)
}

## - Load outlier position lists ──────────────────────────────────────────────
lfmm_pairs    <- read_outlier_pairs(lfmm_path)
pcadapt_pairs <- read_outlier_pairs(pcadapt_path)

## - Normalize chromosome names ───────────────────────────────────────────────
snpspos_norm <- snpspos %>% dplyr::mutate(chrom = normalize_chr(as.character(chrom)))
lfmm_pairs    <- lfmm_pairs    %>% dplyr::mutate(chr = normalize_chr(as.character(chr)))
pcadapt_pairs <- pcadapt_pairs %>% dplyr::mutate(chr = normalize_chr(as.character(chr)))

## - Load genes from GFF and make GRanges (USE ONLY 'gene') ───────────────────
gff <- rtracklayer::import(gff_path)
genes_gr <- gff[gff$type == "gene"]

get_gene_id <- function(gr) {
  md   <- mcols(gr)
  keys <- c("ID","gene_id","Name","gene","locus_tag")
  hit  <- keys[keys %in% colnames(md)][1]
  if (!length(hit)) stop("No recognizable gene ID field in GFF.")
  as.character(md[[hit]])
}
mcols(genes_gr)$gene_id <- get_gene_id(genes_gr)

seqlevels(genes_gr) <- normalize_chr(seqlevels(genes_gr))
genes_gr$chrom      <- normalize_chr(as.character(seqnames(genes_gr)))

## - SNPs → Genes overlap (host gene for each SNP) ────────────────────────────
snps_gr <- GRanges(
  seqnames = normalize_chr(as.character(snpspos_norm$chrom)),
  ranges   = IRanges(start = snpspos_norm$pos, end = snpspos_norm$pos),
  snp_id   = snpspos_norm$snps
)

ov <- findOverlaps(snps_gr, genes_gr, select = "all")

snp2gene <- data.frame(
  snps      = mcols(snps_gr)$snp_id[queryHits(ov)],
  host_gene = as.character(mcols(genes_gr)$gene_id[subjectHits(ov)]),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(host_gene = sub("^(?:gene:|mRNA:)", "", host_gene)) %>%
  dplyr::distinct()

## - Significant eQTL pairs and gene pleiotropy ───────────────────────────────
if ("FDR" %in% names(trans_eQTL2)) {
  trans_eQTL2 <- dplyr::filter(trans_eQTL2, FDR <= FDR_cut)
}

eqtl_pairs <- trans_eQTL2 %>%
  dplyr::distinct(snps, eGene = gene) %>%
  dplyr::inner_join(snp2gene, by = "snps") %>%
  dplyr::mutate(host_gene = as.character(host_gene)) %>%
  dplyr::distinct(snps, eGene, host_gene)

gene_pleiotropy <- eqtl_pairs %>%
  dplyr::count(host_gene, name = "n_eGenes") %>%
  dplyr::mutate(pleiotropic_flag = as.integer(n_eGenes > 1))

## - Map outlier SNPs to SNP IDs, then classify genes ─────────────────────────
lfmm_snps <- lfmm_pairs %>%
  dplyr::inner_join(snpspos_norm, by = c("chr"="chrom","pos"="pos")) %>%
  dplyr::pull(snps) %>% unique()

pcadapt_snps <- pcadapt_pairs %>%
  dplyr::inner_join(snpspos_norm, by = c("chr"="chrom","pos"="pos")) %>%
  dplyr::pull(snps) %>% unique()

gene_class <- snp2gene %>%
  dplyr::mutate(LFMM = snps %in% lfmm_snps,
                PCAdapt = snps %in% pcadapt_snps) %>%
  dplyr::group_by(host_gene) %>%
  dplyr::summarize(LFMM = any(LFMM),
                   PCAdapt = any(PCAdapt),
                   .groups = "drop") %>%
  dplyr::mutate(
    method_class = dplyr::case_when(
      LFMM & !PCAdapt ~ "LFMM_only",
      PCAdapt & !LFMM ~ "PCAdapt_only",
      LFMM &  PCAdapt ~ "Shared",
      TRUE            ~ "Non_outlier"
    ),
    method_class = factor(method_class,
                          levels = c("Non_outlier","LFMM_only","PCAdapt_only","Shared")),
    outlier_flag = as.integer(method_class != "Non_outlier")
  )

## - Assemble gene-level table ──────────────────────────
gene_dat <- gene_pleiotropy %>%
  dplyr::left_join(gene_class, by = "host_gene")

## - Sanity prints ────────────────────────────────────────────────────────────
cat("Genes with pleiotropy info:", nrow(gene_dat), "\n")
print(table(method_class = gene_dat$method_class, useNA = "ifany"))
print(table(pleiotropic = gene_dat$pleiotropic_flag, useNA = "ifany"))

## - Fisher’s exact test: outlier vs pleiotropy ───────────────────────────────
cat("\nExact outlier vs non-outlier (genes):\n")
tab_exact <- with(gene_dat, table(outlier_flag, pleiotropic_flag))
print(tab_exact)
if (all(dim(tab_exact) == c(2,2))) {
  fe_exact <- fisher.test(tab_exact)
  cat("Fisher p (two-sided) =", signif(fe_exact$p.value, 4), "\n")
} else {
  cat("Not a 2×2 table; Fisher exact not applicable here.\n")
}

## - Summary by class ─────────────────────────────────────────────────────────
summary_by_class <- gene_dat %>%
  dplyr::group_by(method_class) %>%
  dplyr::summarize(
    n_genes          = dplyr::n(),
    mean_eGenes      = mean(n_eGenes, na.rm = TRUE),
    median_eGenes    = stats::median(n_eGenes, na.rm = TRUE),
    prop_pleiotropic = mean(pleiotropic_flag == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(method_class, levels = c("Non_outlier","LFMM_only","PCAdapt_only","Shared")))

cat("\n— Mean/median #eGenes per class —\n")
print(summary_by_class)

## write to CSV
write.csv(summary_by_class, "~/Downloads/pleiotropy_summary_by_class.csv", row.names = FALSE)


