setwd("../Argo_redo/")

library(tidyverse)
library(caret)
library(clusterSim)
library(FactoMineR)
library(factoextra)
library(pls)

# read in PCA data
pca <- read_table("argo.eigenvec", col_names = FALSE)
eigenval <- scan("argo.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca <- pca %>%
  filter(!ind %in% c("IID","ARG1820.Aligned.sortedByCoord.out.bam","ARG1805.Aligned.sortedByCoord.out.bam","ARG1834.Aligned.sortedByCoord.out.bam","Ames449.Aligned.sortedByCoord.out.bam","Ames695.Aligned.sortedByCoord.out.bam","arg11B-11.Aligned.sortedByCoord.out.bam","btm5-1.Aligned.sortedByCoord.out.bam","btm7B-14.Aligned.sortedByCoord.out.bam"))

population <- c("Coast", "North", "North", "Coast", "North", "Coast", "Coast", "North", "Coast", "Coast","North", "North", "Coast", "Coast", "Coast", "North", "North", "North", "North")

pca <- data.frame(cbind(pca, population))

# remake data.frame
row.names(pca) <- pca$ind
pca <- pca[order(row.names(pca)),]
n <- row.names(pca)
pca$ind <- NULL

# Read in gene expression data
expression_file_name = paste("normalized_readcounts.csv", sep="")

expr <- read.table(expression_file_name,sep=',', row.names = 1, header = T) 

expr <- expr %>%
  dplyr::select(., -c("arg11B.11", "btm5.1", "btm7B.14"))

expr[1,]

# Find mean expression value for each sample
mean_exp <- data.frame(colMeans(expr)) %>%
  filter()

# Change row names for PCA data
row.names(pca) <- row.names(mean_exp)

# Add mean expression values to PCA dataframe
pca <- pca %>%
  mutate(mean_exp = mean_exp[, 1])

# Convert PC1 and PC2 to numeric for plotting purposes.
pca$PC1 <- as.numeric(pca$PC1)
pca$PC2 <- as.numeric(pca$PC2)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:2, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = population)) + geom_point(size = 6)
b <- b + scale_colour_manual(values = c("brown", "aquamarine4", "darkgrey"))
b <- b + coord_equal() + theme_light()
# c <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
c <- b + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme(
    axis.title.x = element_text(size = 14),   # Increase x-axis label size
    axis.title.y = element_text(size = 14),   # Increase y-axis label size
    axis.text.x = element_text(size = 12),    # Increase x-axis tick number size
    axis.text.y = element_text(size = 12)     # Increase y-axis tick number size
  )



ggsave(plot=c, filename="../argo_plots/PLOS_genetics_figures/pca.tiff", width=7, height=5, device = "tiff", unit = "in")
