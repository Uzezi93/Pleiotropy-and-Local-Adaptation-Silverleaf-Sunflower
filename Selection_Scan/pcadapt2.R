library(pcadapt)
library(qvalue)
library(readr)
library(dplyr)
library(stringr)

# ---------- 1) Read PLINK & run pcadapt ----------
bed_prefix <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/argo_snps"
geno <- pcadapt::read.pcadapt(paste0(bed_prefix, ".bed"), type = "bed")

# Scree to choose K
x0 <- pcadapt(geno, K = 18)         # exploratory
plot(x0, option = "screeplot")      # pick elbow (say K=4 or 5)

# Final run (example: K=5) with LD clumping (adjust as you prefer)
res <- pcadapt(geno, K = 2, LD.clumping = list(size = 200, thr = 0.1))
plot(res)                           # QC plots (stat.distribution, manhattan, qqplot, etc.)


Kset <- c(1,2,3,4)
outs <- list(); lambdas <- c()

for (K in Kset) {
  fit <- pcadapt(matrix, K = K, LD.clumping = list(size = 200, thr = 0.1))
  pv  <- fit$pvalues
  # inflation factor
  lambda <- median(qchisq(1 - pv, 1), na.rm=TRUE) / qchisq(0.5,1)
  lambdas[paste0("K",K)] <- lambda
  qv <- p.adjust(pv, "BH")
  outs[[paste0("K",K)]] <- which(qv < 0.05)
}
lambdas                          # want ~1
lengths(outs)                    # #outliers per K
length(intersect(outs$K1, outs$K2))  # stability

subset_core <- Reduce(intersect, outs[c("K1", "K2","K3","K4")])
length(subset_core)



# Multiple testing
qv <- qvalue(res$pvalues)$qvalues
sum(qv < 0.05, na.rm = TRUE) # how many outliers 3497

# ---------- 2) Population labels in the *correct order* ----------
# Read FAM (order of individuals in bed)
fam <- readr::read_table2(
  paste0(bed_prefix, ".fam"),
  col_names = c("FID","IID","PID","MID","SEX","PHENO"),
  col_types = readr::cols(
    FID   = readr::col_character(),
    IID   = readr::col_character(),
    PID   = readr::col_character(),
    MID   = readr::col_character(),
    SEX   = readr::col_integer(),
    PHENO = readr::col_double()
  )
)


# Make a mapping from fam IDs to your coords table
coords <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/arg_native_poplatlong.csv") |>
  filter(pop %in% pop[1:22]) |>
  filter(!pop %in% c("arg11B","btm5","btm7B")) |>
  mutate(integer = if_else(X=="Coast", 1L, if_else(X=="North", 2L, NA_integer_)))

base_ids <- sub("-.*$", "", fam$IID)                  # strip -4, -6, etc
ord <- match(base_ids, coords$pop)
stopifnot(!any(is.na(ord)))                           # must all match

pop_labels <- coords$X[ord]                           # "Coast"/"North" in fam order
pop_int    <- coords$integer[ord]

# Scores plot with correct order
plot(res, option = "scores", pop = pop_labels)
# plot(res, option = "scores", i = 2, j = 3, pop = pop_labels)

# ---------- 3) Coordinates of outlier SNPs ----------
idx_out <- which(qv < 0.05) #1,796

# Read BIM (same order as bed)
bim <- read_tsv(paste0(bed_prefix, ".bim"),
                col_names = c("CHR","SNP","CM","POS","REF","ALT"),
                show_col_types = FALSE)

stopifnot(nrow(bim) == nrow(res$loadings))            # sanity check order/length

pcadapt_out_coords <- bim[idx_out, c("CHR","POS","SNP","REF","ALT")]
write.csv(pcadapt_out_coords, "pcadapt_outliers_coords.csv", row.names = FALSE)

# Optional: if you prefer variant IDs like chr:pos:ref:alt for downstream
bim_id <- bim |>
  mutate(ID = paste(CHR, POS, REF, ALT, sep=":")) |>
  dplyr::select(ID)
write.table(bim_id$ID[idx_out], "pcadapt_outlier_ids.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#--------------------------------------------------------------------
pcadapt_pos <- bim_id$ID[idx_out] # 1797

#  Compare LFMM and PCAdapt candidates ****
common_loci <- intersect(pcadapt_pos, lfmm_pos) # 435 loci

write.csv(common_loci, "common_outliers.csv", row.names = FALSE)

# Save common SNPs identified from both methods
# create a dataframe for all identified lfmm regions
common_regions <- data.frame(common_loci)

# Split dataframe by ":" and "_" seperators
common_regions <- str_split_fixed(common_regions$common_loci, "_", 2)
common_regions <- str_split_fixed(common_regions[,1], ":", 2)
# common_regions <- str_split_fixed(common_regions[,2], ":", 2)
common_regions1 <- data.frame(common_regions)

# Extract all identified chromosomes
chrom <- common_regions1$X1

common_regions2 <- data.frame(do.call('rbind', strsplit(as.character(common_regions1$X2),'_',fixed=TRUE)))

# Add all separated components into a single dataframe. 
common_region_fin <- cbind(chrom, common_regions2)

# Rename columns
colnames(common_region_fin) <- c("chrm", "pos")

# keep only the numeric position before the first ":"
common_region_fin$pos_clean <- sub(":.*", "", common_region_fin$pos)

# overwrite pos with cleaned values if you prefer
common_region_fin$pos <- common_region_fin$pos_clean
common_region_fin$pos_clean <- NULL


# Save positions rowise. 
write.table(common_region_fin, 
            file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_pcadapt_common_pos.txt", 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")

#--------------------------------------------------------------------
# create a dataframe for all identified pcadapt regions
pcadapt_regions <- data.frame(pcadapt_pos)

# load stringr library
library(stringr)
# Split dataframe by ":" and "_" seperators
pcadapt_regions <- str_split_fixed(pcadapt_regions$pcadapt_pos, "_", 2)
pcadapt_regions <- str_split_fixed(pcadapt_regions[,1], ":", 2)
pcadapt_regions1 <- data.frame(pcadapt_regions)

# Extract all identified chromosomes
chrom <- pcadapt_regions1$X1

pcadapt_regions2 <- data.frame(do.call('rbind', strsplit(as.character(pcadapt_regions1$X2),'_',fixed=TRUE)))

# Add all separated components into a single dataframe. 
pcadapt_region_fin <- cbind(chrom, pcadapt_regions2)

# Rename columns
colnames(pcadapt_region_fin) <- c("chrm", "pos")

# keep only the numeric position before the first ":"
pcadapt_region_fin$pos_clean <- sub(":.*", "", pcadapt_region_fin$pos)

# overwrite pos with cleaned values if you prefer
pcadapt_region_fin$pos <- pcadapt_region_fin$pos_clean
pcadapt_region_fin$pos_clean <- NULL


# Save positions rowise. 
write.table(pcadapt_region_fin, 
            file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/pcadapt_pos.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

#--------------------------------------------------------------------
# Fisher's exact test
# Find common loci
common_loci <- intersect(lfmm_pos, pcadapt_pos) # 435 loci. 
length_common <- length(common_loci)

# Count loci
length_lfmm <- length(lfmm_loci)
length_PCADapt <- length(PCAdapt_loci)
length_not_common <- length_lfmm + length_PCADapt - length_common

# Create contingency table
#        | in common | not in common |
# method1|    10     |               |
# method2|           |               |

# Contingency table
contingency_table <- matrix(c(length_common, length_lfmm - length_common, 
                              length_PCADapt - length_common, length_not_common),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("Method1 (lfmm)", "Method2 (PCAdapt)"),
                                            c("Common", "Not Common")))

# Apply Fisher's exact test
fisher_result <- fisher.test(contingency_table)

# Print results
print(fisher_result) # p-value < 2.2e-16

#-----------------------------------------------------------------------------------------------
# Get pcadapt p-values
pval <- res$pvalues

# Add p-values to your bim table
bim_pval_pcadapt <- bim %>%
  dplyr::select(CHR, POS, SNP, REF, ALT) %>%
  dplyr::mutate(PVAL = pval)

# check
head(bim_pval_pcadapt)

pcadapt_gene_coord <- read.table("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/pcadapt_chrom_pos_gene.txt", sep = "\t")
names(pcadapt_gene_coord) <-  c("CHR", "POS", "gene")

# ensure POS types match (both numeric or integer)
bim_pval_pcadapt <- bim_pval_pcadapt %>%
  mutate(POS = as.integer(POS))   # make sure same type as pcadapt_gene_coord

# join by CHR and POS
merged_df <- bim_pval_pcadapt %>%
  left_join(pcadapt_gene_coord, by = c("CHR", "POS"))

merged_df <- merged_df %>%
  separate(gene, into = c("gene1", "gene2"), sep = ",", fill = "right", remove = TRUE)

# check results
head(merged_df, 10)


# Save positions rowise. 
write.table(merged_df, 
            file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/pcadapt_pvalues.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")



