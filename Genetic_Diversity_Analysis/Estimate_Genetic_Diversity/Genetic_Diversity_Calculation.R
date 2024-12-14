setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/")

# Load required libraries
library(tidyverse)

# Load genetic diversity estimates for all samples calculated with angsd
all_pop_theta <- read.table(file = "all_pop.thetas.idx.pestPG", sep = "\t", header = F)
colnames(all_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

# Estimating population diversity stats
all_pop_theta2 <- all_pop_theta %>%
  dplyr::select(., tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  na.omit()
mean(all_pop_theta2$Watterson_theta) # 0.2225602

all_pop_theta3 <- all_pop_theta %>%
  dplyr::select(., tP, numSites) %>%
  mutate(theta_pie = tP/numSites) %>%
  na.omit()
mean(all_pop_theta3$theta_pie) # 0.2019096

mean(all_pop_theta$tajD) # -0.09140723

all_pop_FST <- read.table(file="all_pop_theta_stat",sep="\t",header=F)
colnames(all_pop_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")
mean(all_pop_FST$weighted_FST) # 0.2133738
