setwd("../argo_plots/")

library(tidyverse)
library(data.table)
library(lubridate)

dxy <- read.table(file="pixy_dxy.txt",sep="\t",header=T);

dxy <- dxy[order(dxy$chromosome),] %>%
  na.omit()

eQTL_pos <- read.table(file="../network_analysis/eQTL_position2.txt",sep="\t",header=F);

# sort position dataframe
eQTL_pos <- with(eQTL_pos, eQTL_pos[order(V1, V2) , ]) %>%
  na.omit

# Ensure chromosome names are consistent in both dataframes
eQTL_pos$V1 <- as.character(eQTL_pos$V1)
dxy$chromosome <- as.character(dxy$chromosome)

# Sort eGenes_pos dataframe by chromosome and position
eQTL_pos <- eQTL_pos %>%
  arrange(V1, V2) %>%
  na.omit()

# Inner join based on matching chromosome names
eQTL_dxy <- inner_join(eQTL_pos, dxy, by = c("V1" = "chromosome"), multiple = "all")

# Now filter for positions where pos is between window_pos_1 and window_pos_2
eQTL_dxy_filtered <- eQTL_dxy %>%
  group_by(V1) %>%
  filter(between(V2, window_pos_1, window_pos_2)) %>%
  select(V1, V2, avg_dxy)

names(eQTL_dxy_filtered) <- c("chrom", "pos", "dxy_stat")

write.table(eQTL_dxy_filtered, file = "eQTL_dxy.txt", quote = F, sep ='\t')

#-------------------------------------------------------------------------
# Extract dxy values for eGenes

eGenes_pos <- read.table(file="../network_analysis/eGenes_position3.txt",sep="\t",header=F);

# sort position dataframe
eGenes_pos <- with(eGenes_pos, eGenes_pos[order(V1, V2) , ]) %>%
  na.omit

# Ensure chromosome names are consistent in both dataframes
eGenes_pos$V1 <- as.character(eGenes_pos$V1)
dxy$chromosome <- as.character(dxy$chromosome)

# Sort eGenes_pos dataframe by chromosome and position
eGenes_pos <- eGenes_pos %>%
  arrange(V1, V2) %>%
  na.omit()

# Inner join based on matching chromosome names
eGenes_dxy <- inner_join(eGenes_pos, dxy, by = c("V1" = "chromosome"), multiple = "all")

# Now filter for positions where pos is between window_pos_1 and window_pos_2
eGenes_dxy_filtered <- eGenes_dxy %>%
  group_by(V1) %>%
  filter(between(V2, window_pos_1, window_pos_2)) %>%
  select(V1, V2, avg_dxy)

names(eGenes_dxy_filtered) <- c("chrom", "pos", "dxy_stat")


write.table(eGenes_dxy_filtered, file = "eGenes_dxy.txt", quote = F, sep ='\t')

#-------------------------------------------------------------------------
# Extract LFMM dxy values
LFMM_pos <- read.table(file="/project/pi_brook_moyers_umb_edu/Uzezi_argo/angsd/FST/lfmm_pos2.txt",sep="\t",header=F) %>%
  separate(V1, into = c("chr", "pos"))

# sort position dataframe
# LFMM_pos <- with(LFMM_pos, LFMM_pos[order(chr, pos) , ]) %>%
#  na.omit

# Ensure chromosome names are consistent in both dataframes
LFMM_pos$chr <- as.character(LFMM_pos$chr)
dxy$chromosome <- as.character(dxy$chromosome)

# Sort LFMM_pos dataframe by chromosome and position
LFMM_pos <- LFMM_pos %>%
  arrange(chr, pos) %>%
  na.omit()

# Inner join based on matching chromosome names
LFMM_dxy <- inner_join(LFMM_pos, dxy, by = c("chr" = "chromosome"), multiple = "all")

# Now filter for positions where pos is between window_pos_1 and window_pos_2
LFMM_dxy_filtered <- LFMM_dxy %>%
  group_by(chr) %>%
  filter(between(pos, window_pos_1, window_pos_2)) %>%
  select(chr, pos, avg_dxy)

names(LFMM_dxy_filtered) <- c("chrom", "pos", "dxy_stat")

write.table(LFMM_dxy_filtered, file = "LFMM_dxy.txt", quote = F, sep ='\t')

#-------------------------------------------------------------------------
# Extract PCAdapt DXY values

PCAdapt_pos <- read.table(file="pcadapt_pos2.txt",sep="\t",header=F) %>%
  separate(V1, into = c("chr", "pos"))

# sort position dataframe
PCAdapt_pos <- with(PCAdapt_pos, PCAdapt_pos[order(chr, pos) , ]) %>%
  na.omit

# Ensure chromosome names are consistent in both dataframes
PCAdapt_pos$chr <- as.character(PCAdapt_pos$chr)
dxy$chromosome <- as.character(dxy$chromosome)

# Sort eGenes_pos dataframe by chromosome and position
PCAdapt_pos <- PCAdapt_pos %>%
  arrange(chr, pos) %>%
  na.omit()

# Inner join based on matching chromosome names
PCAdapt_dxy <- inner_join(PCAdapt_pos, dxy, by = c("chr" = "chromosome"), multiple = "all")

# Now filter for positions where pos is between window_pos_1 and window_pos_2
PCAdapt_dxy_filtered <- PCAdapt_dxy %>%
  group_by(chr) %>%
  filter(between(pos, window_pos_1, window_pos_2)) %>%
  select(chr, pos, avg_dxy)

names(PCAdapt_dxy_filtered) <- c("chrom", "pos", "dxy_stat")

write.table(PCAdapt_dxy_filtered, file = "PCAdapt_dxy.txt", quote = F, sep ='\t')

#------------------------------------------------------------------------
# Get control DXY positions

control_pos <- read.table(file="../../angsd/FST/control_pos3.txt",sep="\t",header=F);

# sort position dataframe
control_pos <- with(control_pos, control_pos[order(V1, V2) , ]) %>%
  na.omit

# Ensure chromosome names are consistent in both dataframes
control_pos$V1 <- as.character(control_pos$V1)
dxy$chromosome <- as.character(dxy$chromosome)

# Sort eGenes_pos dataframe by chromosome and position
control_pos <- control_pos %>%
  arrange(V1, V2) %>%
  na.omit()


# Step 1: Filter dxy to only contain chromosomes in control_pos
chromosomes_of_interest <- unique(control_pos$V1)

# Filter dxy for matching chromosomes
filtered_dxy <- dxy %>%
  filter(chromosome %in% chromosomes_of_interest)

# Step 2: Perform the inner join after filtering to reduce size
# Joining based on chromosome
control_dxy <- inner_join(control_pos, filtered_dxy, by = c("V1" = "chromosome"))

control_dxy_filtered <- control_dxy[control_dxy$V2 >= control_dxy$window_pos_1 & control_dxy$V2 <= control_dxy$window_pos_2, c("V1", "V2", "avg_dxy")]

names(control_dxy_filtered) <- c("chrom", "pos", "dxy_stat")

write.table(control_dxy_filtered, file = "control_dxy", quote = F, sep ='\t')



