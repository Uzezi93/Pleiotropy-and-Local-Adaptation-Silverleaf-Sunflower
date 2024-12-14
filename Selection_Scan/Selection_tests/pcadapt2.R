library(pcadapt)
library(adegenet)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
# BiocManager::install("qvalue")
library(qvalue)
library(tidyverse)

argo <- pcadapt::read.pcadapt("../network_analysis/plink/data_keep2.bed", type = "bed")

# Choose the number of k
x <- pcadapt::pcadapt(input = argo, K = 18)

summary(x)

plot(x, option = "screeplot")

plot(x , option = "manhattan")

# samp <- c("btm9", "btm7B", "btm5", "btm34", "arg2B", "btm10", "arg11B", "arg14B", "btm17", "btm13", "btm30", "btm19", "btm20", "btm26", "btm31", "btm21", "btm22", "arg6B", "btm27", "arg4B", "btm25", "btm32")   

pop.list <- c("Coastal", "North Inland", "North Inland", "Coastal", "North Inland", "Coastal", "Coastal", "North Inland", "Coastal", "Coastal", "North Inland", "North Inland", "Coastal", "Coastal", "Coastal", "North Inland", "North Inland", "North Inland", "North Inland") 

# south <- data.frame(cbind(samp, pop.list)) %>%
#  filter(!pop.list %in% "South Inland")

# Populations for score plot
# 1.	Coastal (samples arg6B-1, btm9-4, btm10-5, btm13-4, btm17-4, btm19-1, btm 20-8, btm21-4, btm22-8)

# 2. North Inland (samples arg14B-7, arg2B-4, arg4B-8, btm25-2, btm26-4, btm27-3, btm30-6, btm31-6, btm32-3, and btm34-6)

# 3. South Inland (samples arg11B-11, btm5-1, btm7B-14)

pop.int <- c(1,2,2,1,2,1,1,2,1,1,2,2,1,1,1,2,2,2,2)

# Score plot
# plot(x, option = "scores", pop = pop.int)

plot(x, option = "scores", pop = pop.list)

# Choose K > 2
plot(x, option = "scores", i = 2, j = 3,  pop = pop.list)

# Controlling for LD
matrix <- argo
res <- pcadapt(matrix, K = 18)
plot(res, option = "screeplot")

# Catell's rule indicates K=4.
res <- pcadapt(matrix, K = 4)
plot(res)

# To evaluate if LD might be an issue for your dataset, we recommend to display the loadings (contributions of each SNP to the PC) and to evaluate if the loadings are clustered in a single or several genomic regions.
par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

# snp thining
res <- pcadapt(matrix, K = 18, LD.clumping = list(size = 200, thr = 0.1))
plot(res, option = "screeplot")

# After SNP thinning, we choose K=4.
res <- pcadapt(matrix, K = 4, LD.clumping = list(size = 200, thr = 0.1))
par(mfrow = c(2, 2))
for (i in 1:3)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

plot(res)

par(mfrow=c(1,1))

plot(res, option = "stat.distribution", K = 4)

plot(res, option = "manhattan", K =4)

hist(res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# Check for outliers
plot(res, option = "qqplot", K = 4)

# q-value correction
padj <- qvalue(res$pvalues)

# convert adjusted p values to q values
qvalues <- padj$qvalues

length(which(qvalues < 0.05)) ## how many SNPs have an FDR < 5%?
# 1747 SNPs were detected by PCAdapt.


# To identify outlier SNPs, load genomic data with PLINK. 

#________________________________________________________________________
# Load genomic data in plink format using adegenet package
a <- read.PLINK("../network_analysis/filtered_vcf_plink.raw", parallel = FALSE)
dim(a)

# Check object structure
str(a)

# Check object class
class(a)

# View object head
head(a)

# Remove all ARGS genotype
a <- a[a$pop!="ARG1820.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="ARG1805.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="ARG1834.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="Ames449.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="Ames695.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="arg11B-11.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="btm5-1.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="btm7B-14.Aligned.sortedByCoord.out.bam"]

a1 <- as.matrix(a)

# Which SNPs have an FDR value < 5%
outlier_SNPs <- a1[, which(qvalues < 0.05)] 
PCAdapt_loci <- colnames(outlier_SNPs) #4968 loci identified by PCAdapt

# Save to .txt file for selection analysis
cat(PCAdapt_loci, file = "PCAdapt.txt")


#  Compare LFMM, RDA and PCAdapt candidates ****
c <- c(PCAdapt_loci, lfmm_loci) ## No common loci identified by all selection tests. 

common_loci <- intersect(lfmm_loci, PCAdapt_loci) # 10 loci. 


# create a dataframe for all identified lfmm regions
pcadapt_regions <- data.frame(PCAdapt_loci)

# load stringr library
library(stringr)
# Split dataframe by ":" and "_" seperators
pcadapt_regions <- str_split_fixed(pcadapt_regions$PCAdapt_loci, "_", 2)
pcadapt_regions <- str_split_fixed(pcadapt_regions[,1], ":", 2)
pcadapt_regions1 <- data.frame(pcadapt_regions)

# Extract all identified chromosomes
chrom <- pcadapt_regions1$X1

pcadapt_regions2 <- data.frame(do.call('rbind', strsplit(as.character(pcadapt_regions1$X2),'_',fixed=TRUE)))

# Add all separated components into a single dataframe. 
pcadapt_region_fin <- cbind(chrom, pcadapt_regions2)


# Rename columns
colnames(pcadapt_region_fin) <- c("chrm", "pos")

# Save positions rowise. 
write.table(pcadapt_region_fin, 
            file = "../network_analysis/pcadapt_pos.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)


# Fisher's exact test
# Find common loci
common_loci <- intersect(lfmm_loci, PCAdapt_loci) # 10 loci. 
length_common <- length(common_loci)

# Count loci
length_lfmm <- length(lfmm_loci)
length_PCADapt <- length(PCAdapt_loci)
length_not_common <- length_lfmm + length_PCADapt - 2 * length_common

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
print(fisher_result)

#-------------------------------------------------------------------------------
# Get PCAdapt loci qvalues
pcadapt_pos <- data.frame(a@loc.names)

# Select odd-numbered rows
pcadapt_pos <- pcadapt_pos[seq(1, nrow(pcadapt_pos), by = 2), ]
pcadapt_pvalues <- data.frame(padj$pvalues)
pcadapt_pvalues$positions <- pcadapt_pos
names(pcadapt_pvalues) <- c("pvalues", "Position")

pcadapt_pvalues <- pcadapt_pvalues 

# Split the 'positions' column by ":" and "_" into multiple parts
split_positions <- str_split_fixed(pcadapt_pvalues$Position, "[:_]", 3)

# Convert the split results to a data frame and rename the columns
split_positions <- data.frame(
  Chromosome = split_positions[, 1],
  Position = split_positions[, 2],
  Genotype = split_positions[, 3]
)

# Combine the split columns with the original 'PC1' column
pcadapt_pvalues_split <- cbind(pcadapt_pvalues[, "pvalues", drop = FALSE], split_positions)

# View the resulting data frame
head(pcadapt_pvalues_split)

all_genes_pos <- read.table("../network_analysis/all_gene_pos2.txt")
names(all_genes_pos) <- c("Chromosome", "Position", "Gene")

# Merge by two columns: "Chromosome" and "Position"
merged_df_pcadapt <- merge(pcadapt_pvalues_split, all_genes_pos, by = c("Chromosome", "Position")) %>%
  select(pvalues, Gene)

# Save positions rowise. 
write.table(merged_df_pcadapt, 
            file = "../network_analysis/pcadapt_pvalues.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)


