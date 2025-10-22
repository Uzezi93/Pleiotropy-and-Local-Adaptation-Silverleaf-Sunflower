library(tidyverse)

# Run angsd for north and coast populations to get pestg files
north_pop_theta <- read.table(file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/genetic_diversity/north_pop.sfs.saf.idx.thetas.idx.pestPG", sep = "\t", header = F)

coast_pop_theta <- read.table(file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/genetic_diversity/coast_pop.sfs.saf.idx.thetas.idx.pestPG", sep = "\t", header = F)

# Rename columns
colnames(north_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")
colnames(coast_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

north_pop <- north_pop_theta %>%
  dplyr::select(., tW, numSites, tP, tajD) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  mutate(Theta_pie = tP/numSites) %>%
  na.omit()
tajD <- mean(north_pop$tajD) # 0.30
Watterson_theta <- mean(north_pop$Watterson_theta) # 0.215
Theta_pie <- mean(north_pop$Theta_pie) # 0.250

coast_pop <- coast_pop_theta %>%
  dplyr::select(., tW, numSites, tP, tajD) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  mutate(Theta_pie = tP/numSites) %>%
  na.omit()
tajD <- mean(coast_pop$tajD) # 0.22
Watterson_theta <- mean(coast_pop$Watterson_theta) # 0.247
Theta_pie <- mean(coast_pop$Theta_pie) # 0.270

# Load genetic diversity estimates for all samples calculated with angsd
all_pop_theta <- read.table(file = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/genetic_diversity/all/all_pop.sfs.saf.idx.thetas.idx.pestPG", sep = "\t", header = F)
colnames(all_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

# Estimating population diversity stats
all_pop_theta2 <- all_pop_theta %>%
  dplyr::select(., tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  na.omit()
mean(all_pop_theta2$Watterson_theta) # 0.2227

all_pop_theta3 <- all_pop_theta %>%
  dplyr::select(., tP, numSites) %>%
  mutate(theta_pie = tP/numSites) %>%
  na.omit()
mean(all_pop_theta3$theta_pie) # 0.2682606

mean(all_pop_theta$tajD) # 0.2275078

all_pop_FST <- read.table(file="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/FST/all_pop_fst_stat",sep="\t",header=F)
colnames(all_pop_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")
mean(all_pop_FST$weighted_FST) # 0.2809224

