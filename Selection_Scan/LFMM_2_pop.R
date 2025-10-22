# LFMM redo without southern populations 
setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/")

# install.packages("sp", repos="http://R-Forge.R-project.org")
# install.packages("sdm")
# install.packages("mapview")
# install.packages("dismo")
# install.packages("rgdal")
# remotes::install_github("wmgeolab/rgeoboundaries")
# install.packages("rasterVis")
# install.packages("adegenet")
# install.packages("LandGenCourse")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("LEA", force = TRUE)

set.seed(1000)

library(remotes)
library(raster)
library(sp)
library(dplyr)
library(tidyr)
library(sdm)
library(mapview)
library(dismo)
library(sf)
library(rgeoboundaries)
library(magrittr)
library(rasterVis)
library(adegenet)
# library(LandGenCourse)
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(KernSmooth)
library(LEA)
library(missMDA)
library(ggplot2)
library(cowplot)
library(pegas)


# Downloading the climate data
r <- getData("worldclim", var="bio", res=5)
r <- r[[c(1:19)]]

# Downloading the boundary of United States
USA <- geoboundaries("United States")

# Extracting temperature and precipitation data of United States
r <- mask(r, as_Spatial(USA))

coords <- read.csv("arg_native_poplatlong.csv")

# Sample coordinates (latitude)
coord2 <- coords %>%
  filter(pop %in% pop[1:22]) %>%
  filter(!pop %in% c("arg11B", "btm5", "btm7B"))

lat <- coord2$lat

#lats <- c(27.87247, 27.25746, 27.14243, 29.64619, 29.27076, 27.87247, 27.25746, 28.72903, 28.10970, 27.91720, 28.69902, 27.61768, 28.26540, 29.19997, 27.67927, 27.75700, 28.04006, 28.44285, 28.61525, 27.85657, 29.35995, 28.44285)

# Sample coordinates (longitude)
lons <- coord2$long

# lons <- c(-97.19174, -97.89252, -98.14900, -97.69089, -97.66214, -97.13564, -97.80295, -97.14055, -97.29214, -97.04207, -97.38758, -97.21676, -97.17168, -97.30904, -97.49628, -97.11933, -97.08013, -97.02723, -97.32674, -97.30904, -97.75291, -97.32563)

# Samples corresponding to arrangement on genomic dataset
sample <- coord2$pop

# sample <- c("btm9", "btm7B", "btm5", "btm34", "arg2B", "btm10", "arg11B", "arg14B", "btm17", "btm13", "btm30", "btm19", "btm20", "btm26", "btm31", "btm21", "btm22", "arg6B", "btm27", "arg4B", "btm25", "btm32")   

# Remove southern samples
all_samples <- data.frame(lat, lons, sample) 

# Creating coordinates
longitude <- as.numeric(all_samples$lons)
latitude <- as.numeric(all_samples$lat)
coords <- data.frame(x = longitude, y = latitude) 

# Converting coordinates to spatial points
points <- SpatialPoints(coords, proj4string = r@crs)

# Crop the raster to the extent of the points with a small buffer
ext <- extent(points) + c(-1, 1, -1, 1)  # Adjust the buffer as needed
r_cropped <- crop(r, ext)

# Extracting values at the points
values <- extract(r_cropped, points)

# Create a data frame for spatial coordinates and bioclim variables
df <- cbind.data.frame(coordinates(points), values)

# Add sample names to data frame
df <- cbind(all_samples, df)

# View data frame
print(df)

# Define color groups for the samples

# Coastal samples (arg6B-1, btm9-4, btm10-5, btm13-4, btm17-4, btm19-1, btm 20-8, btm21-4, btm22-8)

# 2. North Inland samples (arg14B-7, arg2B-4, arg4B-8, btm25-2, btm26-4, btm27-3, btm30-6, btm31-6, btm32-3, and btm34-6)

sample_colors <- c("Coast" =  "brown", "North" = "aquamarine4")

population <- coord2$X


df$population <- population

# Ensure that sample colors cover the filtered samples
# sample_colors <- sample_colors[names(sample_colors) %in% df$population]


# Convert the cropped raster to a data frame for ggplot2
r_df <- as.data.frame(rasterToPoints(r_cropped))

# Plot the cropped raster and points using ggplot2
annual_temp <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = bio1/10), alpha = 0.5) +
  geom_point(data = df, aes(x = lons, y = lat, color = population), size = 3) +
  scale_fill_gradientn(colours = terrain.colors(10), name = "Temperature") +
  scale_color_manual(values = sample_colors, name = "Population") +
  labs(title = "Annual mean temperature (°C)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

annual_prec <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = bio12/10), alpha = 0.5) +
  geom_point(data = df, aes(x = lons, y = lat, color = population), size = 3) +
  scale_fill_gradientn(colours = terrain.colors(10), name = "Precipation") +
  scale_color_manual(values = sample_colors, name = "Population") +
  labs(title = "Annual mean precipitation (mm)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

bioclim_plot <- plot_grid(annual_temp, annual_prec)
 # ggsave("PLOS_genetics_figures/bioclim_plot.tiff", plot = bioclim_plot, device = "tiff", width = 10, height = 5)


#------------------------------------------------------------------------
# Load genomic data in plink format using adegenet package
a <- read.PLINK("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/argo_snps.raw")

# Check object dimensions
dim(a)

# Check object structure
str(a)

# Check object class
class(a)

# View object head
head(a)

#Environmental data
env <- df 
str(env)
#write.table(env, file = "~/env.txt", quote = F, sep = "\t")

# Remove all ARGS genotype and southern populations 
a <- a[a$pop!="ARG1820.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="ARG1805.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="ARG1834.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="Ames449.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="Ames695.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="arg11B-11.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="btm5-1.Aligned.sortedByCoord.out.bam"]
a <- a[a$pop!="btm7B-14.Aligned.sortedByCoord.out.bam"]

# a <- as.matrix(a)

# Write genotypes in the lfmm format
# write.lfmm(a, "genotypes.lfmm")

# Write genotypes in the geno format
#write.geno("genotypes.lfmm", "genotypes2.geno")

# creation of an environment gradient file: gradient.env.
# The .env file contains a single ecological variable
# for each individual.


# Perform PCA on the genlight object
pca_result <- glPca(a)

# Extract the proportion of variance explained by each principal component
var_explained <- pca_result$eig / sum(pca_result$eig)

# Create a data frame for plotting
scree_data <- data.frame(
  Principal_Component = 1:length(var_explained),
  Variance_Explained = var_explained
)

# Plot the screeplot using ggplot2
p <- ggplot(scree_data, aes(x = Principal_Component, y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(color = "red", size = 0.5) +
  geom_point(color = "red", size = 0.5) +
  labs(x = "Principal Component",
       y = "Variance Explained") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6)
  )

ggsave("argo_plots/argo_revision/genotype_screeplot.tiff", plot = p, device = "tiff")


#------------------------------------------------------------------------
# The function snmf() can be run on the data with missing genotypes as follows. The completion of the genotypic matrix is based on estimated ancestry coefficients and ancestral genotype frequencies.

project.missing = snmf("genotypes.lfmm", K = 1:12,
                      entropy = TRUE,
                      repetitions = 10,
                    project = "new")

#plot cross-entropy criterion for all runs in the snmf project
plot(project.missing, col = "blue", pch = 19, cex = 1.2)

ce_min <- sapply(1:10, function(k) min(cross.entropy(project.missing, K = k)))
plot(1:10, ce_min, pch=19); ce_min

# select the run with the lowest cross-entropy value
 best = which.min(cross.entropy(project.missing, K = 3))
 my.colors <- c("tomato", "lightblue",
              "olivedrab", "gold")
 barchart(project.missing, K = 3, run = best,
         border = NA, space = 0,
        col = my.colors,
      xlab = "Individuals",
      ylab = "Ancestry proportions",
     main = "Ancestry matrix") -> bp
 axis(1, at = 1:length(bp$order),
   labels = bp$order, las=1,
   cex.axis = .4)

# The snmf project data can be used to impute the missing data as follows
 impute(project.missing, "genotypes.lfmm",
     method = 'mode', K = 3, run = best)

# Proportion of correct imputation results
dat.imp <- read.lfmm("genotypes.lfmm")
colnames(dat.imp) <- a@loc.names

#----------------------------------------------------------------------------------

pred <- env[, 4:24]

nb <- estim_ncpPCA(pred,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)

nb$ncp #4

res.comp <- imputePCA(pred, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
imp <- cbind.data.frame(res.comp$completeObs)

write.env(imp, "gradients.env")


#----------------------------------------------------------------------------------

pred.pca <- rda(imp, scale=T)
plot(pred.pca)
summary(pred.pca)$cont
x <- data.frame(summary(pred.pca)$cont)

y <- x %>%
  mutate_if(is.numeric, round, digits = 2)

 write.table(y, file = "BIOCLIM_PCs.txt", sep = "\t", quote = F)

screeplot(pred.pca, main = "Screeplot: Eigenvalues of Argo Predictor Variables")

## correlations between the PC axis and predictors:
# species (column) scores, first up to 17 axes
# sample scores (rows)
si <- vegan::scores(pred.pca, display = "sites", scaling = 1)
round(si[, 1:min(17, ncol(si)), drop = FALSE], 3)

# get site scores, PC1 only, scaling 0
pred.PC1 <- vegan::scores(pred.pca, display = "sites", choices = 1, scaling = 0)

# as a named vector (often handy)
pred.PC1_vec <- setNames(pred.PC1[,1], rownames(pred.PC1))

screeplot(pred.pca, main = "Screeplot of H.argophyllus Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")

# Genetic scree plot
gen.pca <- rda(dat.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

#------------------------------------------------------------------------
# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 2

argo.lfmm <- lfmm_ridge(Y=dat.imp, X=pred.PC1, K=K) ## change K as you see fit

argo.pv <- lfmm_test(Y=dat.imp, X=pred.PC1, lfmm=argo.lfmm, calibrate="gif")
names(argo.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
argo.pv$gif     #1.453028    (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv$pvalue, main="Unadjusted p-values")        
hist(argo.pv$calibrated.pvalue, main="GIF-adjusted p-values")      
#------------------------------------------------------------------------
# How to manually adjust p-values:
# Let's change the GIF and readjust the p-values:
# zscore <- argo.pv$score[,1:17]   # zscores for first predictor, we only have one in our case...
# (gif <- argo.pv$gif)       ## default GIF for this predictor

# Adjusted to 0.5
# new.gif <- 0.5               ## choose your new GIF

# adj.pv1 <- pchisq(zscore^2/new.gif, df=1, lower = FALSE)

# hist(argo.pv$pvalue[,1], main="Unadjusted p-values")        
# hist(argo.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=3.0)")
# hist(adj.pv1, main="REadjusted p-values (GIF=2.0)")

#------------------------------------------------------------------------
# convert adjusted p values to q values
argo.qv <- qvalue(argo.pv$calibrated.pvalue)$qvalues

length(which(argo.qv < 0.05)) ## how many SNPs have an FDR < 5%? #0

# Using K=1 and default GIF calculated from lfmmm, and an FDR threshold of 0.5, we only detected 151 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.1 <- a[, which(argo.qv < 0.05)] 
#argo.FDR <- a[, which(argo.qv > 0)] 

# Get loci and qvalues for regression
argo.pv1 <- data.frame(argo.pv$pvalue)

loci_and_pval1 <- argo.pv1


#------------------------------------------------------------------------
# For PC2

pred.PC2 <- vegan::scores(pred.pca, display = "sites", choices = 2, scaling = 0)

# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 2

argo.lfmm2 <- lfmm_ridge(Y=dat.imp, X=pred.PC2, K=K) ## change K as you see fit

argo.pv2 <- lfmm_test(Y=dat.imp, X=pred.PC2, lfmm=argo.lfmm2, calibrate="gif")
names(argo.pv2) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
argo.pv2$gif     #1.352139 (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv2$pvalue, main="Unadjusted p-values")        
hist(argo.pv2$calibrated.pvalue, main="GIF-adjusted p-values") 

# convert adjusted p values to q values
argo.qv2 <- qvalue(argo.pv2$calibrated.pvalue)$qvalues

length(which(argo.qv2 < 0.05)) ## how many SNPs have an FDR < 5%? #3706

# Using K=3 and default GIF calculated from lfmmm, and an FDR threshold of 0.05, we only detected 9 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.2 <- a[, which(argo.qv2 < 0.05)] 

# Get loci and qvalues for regression
argo.pv2 <- data.frame(argo.pv2$pvalue)

loci_and_pval2 <- argo.pv2 

#------------------------------------------------------------------------

# For PC3

pred.PC3 <- vegan::scores(pred.pca, display = "sites", choices = 3, scaling = 0)


# Running LFMM (Univariate GEA)
# Run lfmm at K=1
K <- 2

argo.lfmm3 <- lfmm_ridge(Y=dat.imp, X=pred.PC3, K=K) ## change K as you see fit

argo.pv3 <- lfmm_test(Y=dat.imp, X=pred.PC3, lfmm=argo.lfmm3, calibrate="gif")
names(argo.pv3) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
# argo.pv3$gif     #1.138296    (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv3$pvalue, main="Unadjusted p-values")        
hist(argo.pv3$calibrated.pvalue, main="GIF-adjusted p-values") 

# convert adjusted p values to q values
argo.qv3 <- qvalue(argo.pv3$calibrated.pvalue)$qvalues

length(which(argo.qv3 < 0.05)) ## how many SNPs have an FDR < 5%? #2358

# Using K=3 and default GIF calculated from lfmmms, and an FDR threshold of 0.5, we only detected 78 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.3 <- a[, which(argo.qv3 < 0.05)] 

# Get loci and pvalues for regression
argo.pv3 <- data.frame(argo.pv3$pvalue)

loci_and_pval3 <- argo.pv3 

#----------------------------------------------------------------------------------

# For PC4

pred.PC4 <- vegan::scores(pred.pca, display = "sites", choices = 4, scaling = 0)


# Running LFMM (Univariate GEA)
# Run lfmm at K=1
K <- 2

argo.lfmm4 <- lfmm_ridge(Y=dat.imp, X=pred.PC4, K=K) ## change K as you see fit

argo.pv4 <- lfmm_test(Y=dat.imp, X=pred.PC4, lfmm=argo.lfmm4, calibrate="gif")
names(argo.pv4) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
# argo.pv3$gif     #1.138296    (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv4$pvalue, main="Unadjusted p-values")        
hist(argo.pv4$calibrated.pvalue, main="GIF-adjusted p-values") 

# convert adjusted p values to q values
argo.qv4 <- qvalue(argo.pv4$calibrated.pvalue)$qvalues

length(which(argo.qv4 < 0.05)) ## how many SNPs have an FDR < 5%? #14

# Using K=3 and default GIF calculated from lfmmms, and an FDR threshold of 0.5, we only detected 78 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.4 <- a[, which(argo.qv4 < 0.05)] 

# Get loci and pvalues for regression
argo.pv4 <- data.frame(argo.pv4$pvalue)

loci_and_pval4 <- argo.pv4 


#----------------------------------------------------------------------------------
lfmm_loci1 <- which(argo.qv < 0.05) #0
lfmm_loci2 <- which(argo.qv2 < 0.05) #3706
lfmm_loci3 <- which(argo.qv3 < 0.05) #2358
lfmm_loci4 <- which(argo.qv4 < 0.05) #14
# lfmm_loci5 <- argo.FDR.5@loc.names #99

lfmm_loci <- c(lfmm_loci1, lfmm_loci2, lfmm_loci3, lfmm_loci4) #6078

# ---------- 3) Coordinates of outlier SNPs ----------
idx_out <- lfmm_loci

# Read BIM (same order as bed)
bim <- read_tsv(paste0(bed_prefix, ".bim"),
                col_names = c("CHR","SNP","CM","POS","REF","ALT"),
                show_col_types = FALSE)

lfmm_out_coords <- bim[idx_out, c("CHR","POS","SNP","REF","ALT")]
write.csv(lfmm_out_coords, "lfmm_outliers_coords.csv", row.names = FALSE)

# Optional: if you prefer variant IDs like chr:pos:ref:alt for downstream
bim_id <- bim |>
  mutate(ID = paste(CHR, POS, REF, ALT, sep=":")) |>
  dplyr::select(ID)
write.table(bim_id$ID[idx_out], "lfmm_outlier_ids.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

lfmm_pos <- bim_id$ID[idx_out]


# create a dataframe for all identified lfmm regions
lfmm_regions <- data.frame(lfmm_pos)

# load stringr library
library(stringr)
# Split dataframe by ":" and "_" seperators
lfmm_regions <- str_split_fixed(lfmm_regions$lfmm_pos, "_", 2)
lfmm_regions <- str_split_fixed(lfmm_regions[,1], ":", 2)
lfmm_regions1 <- data.frame(lfmm_regions)

# Extract all identified chromosomes
chrom <- lfmm_regions1$X1

lfmm_regions2 <- data.frame(do.call('rbind', strsplit(as.character(lfmm_regions1$X2),'_',fixed=TRUE)))

# Add all separated components into a single dataframe. 
lfmm_region_fin <- cbind(chrom, lfmm_regions2)

# Rename columns
colnames(lfmm_region_fin) <- c("chrm", "pos")

# keep only the numeric position before the first ":"
lfmm_region_fin$pos_clean <- sub(":.*", "", lfmm_region_fin$pos)

# overwrite pos with cleaned values if you prefer
lfmm_region_fin$pos <- lfmm_region_fin$pos_clean
lfmm_region_fin$pos_clean <- NULL


# Save positions rowise. 
write.table(lfmm_region_fin, 
            file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_pos.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)


# Get LFMM p values
pval1 <- argo.pv1$PC1
pval2 <- argo.pv2$PC2
pval3 <- argo.pv3$PC3
pval4 <- argo.pv4$PC4


# Add p-values to your bim table
# Combine all p-values into one data frame with SNP info
bim_pval_all <- bim %>%
  dplyr::select(CHR, POS, SNP, REF, ALT) %>%
  dplyr::mutate(
    PVAL_PC1 = pval1,
    PVAL_PC2 = pval2,
    PVAL_PC3 = pval3,
    PVAL_PC4 = pval4
  )

# check structure
str(bim_pval_all)
head(bim_pval_all)


bim_pval_all <- bim_pval_all %>%
  mutate(across(starts_with("PVAL_PC"), ~ as.vector(.x)))


lfmm_gene_coord <- read.table("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_chrom_pos_gene.txt", sep = "\t")
names(lfmm_gene_coord) <-  c("CHR", "POS", "gene")

# ensure POS types match (both numeric or integer)
bim_pval_lfmm <- bim_pval_all %>%
  mutate(POS = as.integer(POS))   # make sure same type as pcadapt_gene_coord

# join by CHR and POS
merged_df <- bim_pval_all %>%
  left_join(lfmm_gene_coord, by = c("CHR", "POS"))

merged_df <- merged_df %>%
  separate(gene, into = c("gene1", "gene2"), sep = ",", fill = "right", remove = TRUE)

# check results
head(merged_df, 10)

# Save positions rowise. 
write.table(merged_df, 
            file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/lfmm_pval.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")
