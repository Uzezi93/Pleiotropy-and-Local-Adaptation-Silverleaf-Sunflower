# LFMM redo without southern populations 
setwd("~/Helianthus_argophyllus_project/Argo_redo/")

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

# Sample coordinates (latitude)
lats <- c(27.87247, 27.25746, 27.14243, 29.64619, 29.27076, 27.87247, 27.25746, 28.72903, 28.10970, 27.91720, 28.69902, 27.61768, 28.26540, 29.19997, 27.67927, 27.75700, 28.04006, 28.44285, 28.61525, 27.85657, 29.35995, 28.44285)

# Sample coordinates (longitude)
lons <- c(-97.19174, -97.89252, -98.14900, -97.69089, -97.66214, -97.13564, -97.80295, -97.14055, -97.29214, -97.04207, -97.38758, -97.21676, -97.17168, -97.30904, -97.49628, -97.11933, -97.08013, -97.02723, -97.32674, -97.30904, -97.75291, -97.32563)

# Samples corresponding to arrangement on genomic dataset
sample <- c("btm9", "btm7B", "btm5", "btm34", "arg2B", "btm10", "arg11B", "arg14B", "btm17", "btm13", "btm30", "btm19", "btm20", "btm26", "btm31", "btm21", "btm22", "arg6B", "btm27", "arg4B", "btm25", "btm32")   

# Remove southern samples
all_samples <- data.frame(lats, lons, sample) %>%
  filter(!sample %in% c("arg11B", "btm5", "btm7B"))

# Creating coordinates
longitude <- as.numeric(all_samples$lons)
latitude <- as.numeric(all_samples$lats)
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

# Coastal (samples arg6B-1, btm9-4, btm10-5, btm13-4, btm17-4, btm19-1, btm 20-8, btm21-4, btm22-8)

# 2. North Inland (samples arg14B-7, arg2B-4, arg4B-8, btm25-2, btm26-4, btm27-3, btm30-6, btm31-6, btm32-3, and btm34-6)

sample_colors <- c("Coast" =  "brown", "North" = "aquamarine4")

population <- c("Coast", "North", 
                "North", "Coast", "North", 
                "Coast", "Coast", "North", "Coast", 
                "Coast", "North", "North", "Coast", 
                "Coast", "Coast", "North", "North", 
                "North", "North")

df$population <- population

# Ensure that sample colors cover the filtered samples
# sample_colors <- sample_colors[names(sample_colors) %in% df$population]


# Convert the cropped raster to a data frame for ggplot2
r_df <- as.data.frame(rasterToPoints(r_cropped))

# Plot the cropped raster and points using ggplot2
annual_temp <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = bio1/10), alpha = 0.5) +
  geom_point(data = df, aes(x = lons, y = lats, color = population), size = 3) +
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
  geom_point(data = df, aes(x = lons, y = lats, color = population), size = 3) +
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
a <- read.PLINK("../network_analysis/filtered_vcf_plink.raw")

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
# write.geno("./genotypes.lfmm_imputed.lfmm", "genotypes2.geno")

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
# p <- ggplot(scree_data, aes(x = Principal_Component, y = Variance_Explained)) +
#  geom_bar(stat = "identity", fill = "steelblue") +
 # geom_line(color = "red", size = 0.5) +
#  geom_point(color = "red", size = 0.5) +
 # labs(x = "Principal Component",
  #     y = "Variance Explained") +
#  theme_minimal() +
 # theme(
  #  axis.text = element_text(size = 6),
   # axis.title = element_text(size = 6)
#  )

# ggsave("PLOS_genetics_figures/genotype_screeplot.tiff", plot = p, device = "tiff")


#------------------------------------------------------------------------
# The function snmf() can be run on the data with missing genotypes as follows. The completion of th genotypic matrix is based on estimated ancestry coefficients and ancestral genotype frequencies.

# project.missing = snmf("genotypes.lfmm", K = 1:12,
 #                      entropy = TRUE,
  #                    repetitions = 10,
   #                  project = "new")

# plot cross-entropy criterion for all runs in the snmf project
# plot(project.missing, col = "blue", pch = 19, cex = 1.2)

# select the run with the lowest cross-entropy value
# best = which.min(cross.entropy(project.missing, K = 3))
# my.colors <- c("tomato", "lightblue",
#              "olivedrab", "gold")
# barchart(project.missing, K = 3, run = best,
#         border = NA, space = 0,
 #       col = my.colors,
  #     xlab = "Individuals",
   #   ylab = "Ancestry proportions",
    # main = "Ancestry matrix") -> bp
# axis(1, at = 1:length(bp$order),
  # labels = bp$order, las=1,
  # cex.axis = .4)

# The snmf project data can be used to impute the missing data as follows
# impute(project.missing, "genotypes.lfmm",
#     method = 'mode', K = 3, run = best)

# Proportion of correct imputation results
dat.imp <- read.lfmm("../Argo_redo/genotypes.lfmm_imputed.lfmm")
colnames(dat.imp) <- a@loc.names

#----------------------------------------------------------------------------------

pred <- env[, 4:22]

nb <- estim_ncpPCA(pred,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)

nb$ncp #4

res.comp <- imputePCA(pred, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set
imp <- cbind.data.frame(res.comp$completeObs)

# write.env(imp, "gradients.env")


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
round(scores(pred.pca, choices=1:17, display=c("species", scaling=0), digits=3))

pred.PC1 <- scores(pred.pca, choices=1, display=c("sites", scaling=0))
screeplot(pred.pca, main = "Screeplot of H.argophyllus Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")

# Genetic scree plot
gen.pca <- rda(dat.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

#------------------------------------------------------------------------
# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 4

argo.lfmm <- lfmm_ridge(Y=dat.imp, X=pred.PC1, K=K) ## change K as you see fit

argo.pv <- lfmm_test(Y=dat.imp, X=pred.PC1, lfmm=argo.lfmm, calibrate="gif")
names(argo.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
argo.pv$gif     #1.301588   (Tests seems appropriately calibrated.

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

length(which(argo.qv < 0.05)) ## how many SNPs have an FDR < 5%? #32

# Using K=1 and default GIF calculated from lfmmm, and an FDR threshold of 0.5, we only detected 151 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.1 <- a[, which(argo.qv < 0.05)] 
#argo.FDR <- a[, which(argo.qv > 0)] 

argo.pv1 <- data.frame(argo.pv$pvalue)

loci_and_pval <- argo.pv1 

#------------------------------------------------------------------------
# For PC2

pred.PC2 <- scores(pred.pca, choices=2, display=c("sites", scaling=0))

# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 4

argo.lfmm2 <- lfmm_ridge(Y=dat.imp, X=pred.PC2, K=K) ## change K as you see fit

argo.pv2 <- lfmm_test(Y=dat.imp, X=pred.PC2, lfmm=argo.lfmm, calibrate="gif")
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

length(which(argo.qv2 < 0.05)) ## how many SNPs have an FDR < 5%? #162

# Using K=3 and default GIF calculated from lfmmm, and an FDR threshold of 0.05, we only detected 9 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.2 <- a[, which(argo.qv2 < 0.05)] 

# Get loci and qvalues for regression
argo.pv2 <- data.frame(argo.pv2$pvalue)

loci_and_pval2 <- argo.pv2 

#------------------------------------------------------------------------

# For PC3

pred.PC3 <- scores(pred.pca, choices=3, display=c("sites", scaling=0))

# Running LFMM (Univariate GEA)
# Run lfmm at K=1
K <- 4

argo.lfmm3 <- lfmm_ridge(Y=dat.imp, X=pred.PC3, K=K) ## change K as you see fit

argo.pv3 <- lfmm_test(Y=dat.imp, X=pred.PC3, lfmm=argo.lfmm, calibrate="gif")
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

length(which(argo.qv3 < 0.05)) ## how many SNPs have an FDR < 5%? #1192

# Using K=3 and default GIF calculated from lfmmms, and an FDR threshold of 0.5, we only detected 78 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.3 <- a[, which(argo.qv3 < 0.05)] 

# Get loci and qvalues for regression
argo.qv3 <- data.frame(argo.qv3)

loci_and_qval3 <- argo.qv3 

#------------------------------------------------------------------------
# For PC4

pred.PC4 <- scores(pred.pca, choices=4, display=c("sites", scaling=0))

# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 4

argo.lfmm4 <- lfmm_ridge(Y=dat.imp, X=pred.PC4, K=K) ## change K as you see fit

argo.pv4 <- lfmm_test(Y=dat.imp, X=pred.PC4, lfmm=argo.lfmm, calibrate="gif")
names(argo.pv4) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
argo.pv4$gif     #1.227405  (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv4$pvalue, main="Unadjusted p-values")        
hist(argo.pv4$calibrated.pvalue, main="GIF-adjusted p-values") 

# convert adjusted p values to q values
argo.qv4 <- qvalue(argo.pv4$calibrated.pvalue)$qvalues

length(which(argo.qv4 < 0.05)) ## how many SNPs have an FDR < 5%? #26

# Using K=3 and default GIF calculated from lfmmms, and an FDR threshold of 0.5, we only detected 78 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.4 <- a[, which(argo.qv4 < 0.05)]

# Get loci and qvalues for regression
argo.qv4 <- data.frame(argo.qv4)

loci_and_qval4 <- argo.qv4 

#--------------------------------------------------------------------------------
# For PC5

pred.PC5 <- scores(pred.pca, choices=5, display=c("sites", scaling=0))

# Running LFMM (Univariate GEA)
# Run lfmm at K=4
K <- 4

argo.lfmm5 <- lfmm_ridge(Y=dat.imp, X=pred.PC5, K=K) ## change K as you see fit

argo.pv5 <- lfmm_test(Y=dat.imp, X=pred.PC5, lfmm=argo.lfmm, calibrate="gif")
names(argo.pv5) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

# look at the Genomic Inflation Factor(GIF)
argo.pv5$gif     #1.227405  (Tests seems appropriately calibrated.

# An appropriately calibrated set of tests will have a GIF of around 1. An elevated GIF indicates that the results may be overly liberal in identifying candidate SNPs. 
#If the GIF is less than one, the test may be too conservative.

# How application of the GIF to the p-values impacts the p-value distribution:

hist(argo.pv5$pvalue, main="Unadjusted p-values")        
hist(argo.pv5$calibrated.pvalue, main="GIF-adjusted p-values") 

# convert adjusted p values to q values
argo.qv5 <- qvalue(argo.pv5$calibrated.pvalue)$qvalues

length(which(argo.qv5 < 0.05)) ## how many SNPs have an FDR < 5%? #99

# Using K=3 and default GIF calculated from lfmmms, and an FDR threshold of 0.5, we only detected 78 candidate SNPs under selection in response to our PC1 environmental predictor.

## identify which SNPs these are
argo.FDR.5 <- a[, which(argo.qv4 < 0.05)]

# Get loci and qvalues for regression
argo.qv5 <- data.frame(argo.qv5)

loci_and_qval5 <- argo.qv5 

#----------------------------------------------------------------------------------
lfmm_loci1 <- colnames(argo.FDR.1) #44
lfmm_loci2 <- colnames(argo.FDR.2) #125
lfmm_loci3 <- colnames(argo.FDR.3) #1271
lfmm_loci4 <- colnames(argo.FDR.4) #15
lfmm_loci5 <- colnames(argo.FDR.5) #138

lfmm_loci <- c(lfmm_loci1, lfmm_loci2, lfmm_loci3, lfmm_loci4, lfmm_loci5) #1438

# Save to .txt file for selection analysis
cat(lfmm_loci, file = "lfmm3.txt")

# create a dataframe for all identified lfmm regions
lfmm_regions <- data.frame(lfmm_loci)

# load stringr library
library(stringr)
# Split dataframe by ":" and "_" seperators
lfmm_regions <- str_split_fixed(lfmm_regions$lfmm_loci, ":", 2)
# sepearting by "_"
lfmm_regions1 <- data.frame(lfmm_regions)

# Extract all identified chromosomes
chrom <- lfmm_regions1$X1

lfmm_regions2 <- data.frame(do.call('rbind', strsplit(as.character(lfmm_regions1$X2),'_',fixed=TRUE)))

# Add all separated components into a single dataframe. 
lfmm_region_fin <- cbind(chrom, lfmm_regions2)

# Rename columns
colnames(lfmm_region_fin) <- c("chrm", "pos","SNP")

# save positions for extracting positions from VCF file
cat(lfmm_region_fin$pos, file = "lfmm3_pos.txt")

# Save positions rowise. 
write.table(lfmm_region_fin, 
            file = "../network_analysis/lfmm3_pos.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)



# Get lfmm positions and corresponding qvalues
colnames(loci_and_qval2) <- colnames(loci_and_qval)
colnames(loci_and_qval3) <- colnames(loci_and_qval)
colnames(loci_and_qval4) <- colnames(loci_and_qval)
colnames(loci_and_qval5) <- colnames(loci_and_qval)

lfmm_qvalues <- rbind(loci_and_qval, loci_and_qval2, loci_and_qval3, loci_and_qval4, loci_and_qval5)
lfmm_qvalues$positions <- rownames(lfmm_qvalues)




lfmm_pvalues <- loci_and_pval
lfmm_pvalues$positions <- rownames(lfmm_pvalues)
# Split the 'positions' column by ":" and "_" into multiple parts
split_positions <- str_split_fixed(lfmm_pvalues$positions, "[:_]", 3)

# Convert the split results to a data frame and rename the columns
split_positions <- data.frame(
  Chromosome = split_positions[, 1],
  Position = split_positions[, 2],
  Genotype = split_positions[, 3]
)

# Combine the split columns with the original 'PC1' column
lfmm_pvalues_split <- cbind(lfmm_pvalues[, "PC1", drop = FALSE], split_positions)

# View the resulting data frame
head(lfmm_pvalues_split)

all_genes_pos <- read.table("../network_analysis/all_gene_pos2.txt")
names(all_genes_pos) <- c("Chromosome", "Position", "Gene")

# Merge by two columns: "Chromosome" and "Position"
merged_df <- merge(lfmm_pvalues_split, all_genes_pos, by = c("Chromosome", "Position")) %>%
  select(PC1, Gene)

names(merged_df) <- c("pvalues", "Gene")

# Save positions rowise. 
write.table(merged_df, 
            file = "../network_analysis/lfmm_pvalues.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)














