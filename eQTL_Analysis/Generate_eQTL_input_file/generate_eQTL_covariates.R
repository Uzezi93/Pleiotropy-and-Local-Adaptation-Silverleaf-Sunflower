# Load library
library(tidyverse)

# Get covariate for eQTL analysis (K = 4)
#read in the data for k=4
q <-read.table("~/admixture/str_k_4.qopt")

#read in the population information file located in the bam.filelist
pop<-read.table("~/admixture/samples.txt")

q <- data.frame(cbind(q,pop))
names(q) <- c("K1","K2","K3","K4","individuals")

q <- q %>%
  filter(!individuals %in% c("ARG1820.Aligned.sortedByCoord.out.bam", "ARG1805.Aligned.sortedByCoord.out.bam", "ARG1834.Aligned.sortedByCoord.out.bam", "Ames695.Aligned.sortedByCoord.out.bam", "Ames449.Aligned.sortedByCoord.out.bam"))

covariates <- data.frame(t(q))

# make header sample names
# Extract the last row
new_header <- as.character(covariates[nrow(covariates), ])

# Remove the last row from the dataframe
covariates <- covariates[-nrow(covariates), ]

# Assign the new header
colnames(covariates) <- new_header
covariates$`btm9-4.Aligned.sortedByCoord.out.bam` <- NULL

write.table(covariates, file = "../eQTL_analysis/covariates.txt", sep = "\t", quote = F)
