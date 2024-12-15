setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/")

# Plot CLR for PCadapt loci for north and coast populations
library(stringr)
library(tidyverse)
library(tidyr)
library(data.table)
library(tibble)
library(formattable)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(rstatix)
library(ggsignif)


north_lfmm <-read.table(file="north_lfmm_clr.out",sep="\t",header=T);

coast_lfmm <-read.table(file="coast_lfmm_clr.out",sep="\t",header=T);

north_pcadapt <- read.table(file="north_pcadapt_clr.out",sep="\t",header=T);

coast_pcadapt <- read.table(file="coast_pcadapt_clr.out",sep="\t",header=T);

north_eQTL <- read.table(file="north_eQTL2_clr.out",sep="\t",header=T);

coast_eQTL <- read.table(file="coast_eQTL2_clr.out",sep="\t",header=T);

north_eGenes <- read.table(file="north_eGenes2_clr.out",sep="\t", fill = T, header=T);

coast_eGenes <- read.table(file="coast_eGenes2_clr.out",sep="\t",fill = T, header=T);

north_control <-read.table(file="north_control2_clr.out",sep="\t",header=T);

coast_control <- read.table(file="coast_control2_clr.out",sep="\t",header=T);

coast_shared <- read.table(file = "coast_shared_clr.out", sep = "\t", header = T)

north_shared <- read.table(file = "north_shared_clr.out", sep = "\t", header = T)
#------------------------------------------------------------------------
# Load theta and neutrality stats
north_lfmm_theta <- read.table(file="../../angsd/SFS2/north_lfmm.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

north_PCAdapt_theta <- read.table(file="../../angsd/SFS2/north_pcadapt.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

coast_lfmm_theta <- read.table(file="../../angsd/SFS2/coast_lfmm.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

coast_PCAdapt_theta <- read.table(file="../../angsd/SFS2/coast_pcadapt.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

north_eQTL_theta <- read.table(file="../../angsd/SFS3/north_eQTL.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

coast_eQTL_theta <- read.table(file="../../angsd/SFS3/coast_eQTL.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

north_eGenes_theta <- read.table(file="../../angsd/SFS3/north_eGenes.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

coast_eGenes_theta <- read.table(file="../../angsd/SFS3/coast_eGenes.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

north_control_theta <- read.table(file="../../angsd/SFS3/north_control.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

coast_control_theta <- read.table(file="../../angsd/SFS3/coast_control.sfs.saf.idx.thetas.idx.pestPG",sep="\t",header=F);

all_pop_theta <- read.table(file = "../../angsd/SFS2/all_pop.thetas.idx.pestPG", sep = "\t", header = F)

coast_shared_theta <- read.table(file = "../../angsd/SFS2/coast_shared.thetas.idx.pestPG", sep = "\t", header = F)

north_shared_theta <- read.table(file = "../../angsd/SFS2/north_shared.thetas.idx.pestPG", sep = "\t", header = F)


# Rename columns
colnames(north_lfmm_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_lfmm_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(north_PCAdapt_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_PCAdapt_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_eGenes_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(north_eGenes_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(north_eQTL_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_eQTL_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(north_control_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_control_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(all_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(coast_shared_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

colnames(north_shared_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#------------------------------------------------------------------------
# Load in FST data
Control_FST <- read.table(file="../../angsd/FST/control_fst.txt",sep="\t",header=F);

PCAdapt_FST <- read.table(file="../../angsd/FST/pcadapt_fst.txt",sep="\t",header=F);

LFMM_FST <- read.table(file="../../angsd/FST/lfmm_fst.txt",sep="\t",header=F);

eQTL_FST <- read.table(file="../../angsd/FST/eQTL_fst.txt",sep="\t",header=F);

eGene_FST <- read.table(file="../../angsd/FST/eGenes_fst.txt",sep="\t",header=F);

#Rename columns for FST data
colnames(Control_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

colnames(PCAdapt_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

colnames(LFMM_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

colnames(eQTL_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

colnames(eGene_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

#------------------------------------------------------------------------
# Load DXY Data
Control_dxy <- read.table(file="control_dxy",sep="\t",header=T);

LFMM_dxy <- read.table(file="LFMM_dxy.txt",sep="\t",header=T);

PCAdapt_dxy <- read.table(file="PCAdapt_dxy.txt",sep="\t",header=T);

eQTL_dxy <- read.table(file="eQTL_dxy.txt",sep="\t",header=T);

eGene_dxy <- Control_theta <- read.table(file="eGenes_dxy.txt",sep="\t",header=T);

#------------------------------------------------------------------------
# Plot CLR distributions Northern populations
north_clr_LFMM <- north_lfmm$LR
north_clr_PCAdapt <- north_pcadapt$LR
north_clr_eQTL <- north_eQTL$LR
north_clr_eGene <- north_eGenes$LR
north_clr_Control <- north_control$LR

north_clr_data <- cbind(north_clr_LFMM, north_clr_PCAdapt, north_clr_eQTL, north_clr_eGene, north_clr_Control)
north_clr_data <- as.data.frame(north_clr_data) 
names(north_clr_data) <- c("LFMM", "PCAdapt", "eQTL", "eGene", "Control")

north_clr_data_long <- north_clr_data %>%
  gather(Method, CLR) %>%
  mutate(CLR = log10(CLR + 1))

north_clr_data_long <- north_clr_data_long[is.finite(north_clr_data_long$CLR), ]

median(filter(north_clr_data_long, Method == "Control")$CLR) # 0.01600694


# Define color palette for groups
my_colors <- c("Control" = "darkgrey", "eGene" = "brown", "eQTL" = "brown", "LFMM" = "aquamarine4", "PCAdapt" = "aquamarine4")

p <- north_clr_data_long %>%
  ggplot(aes(x=Method, y=CLR, fill=Method)) +
  geom_violin(adjust = 10, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.25, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.077, color="black", alpha = 0.8) +
  scale_y_continuous(trans = 'log10') +
  xlab(" ") +
  ylab("CLR") +
  ggtitle("North") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_CLR.tiff", plot = p, device = "tiff")

clr_lfmm <- filter(north_clr_data_long, Method == "LFMM")$CLR
clr_pcad <- filter(north_clr_data_long, Method == "PCAdapt")$CLR
clr_control <- filter(north_clr_data_long, Method == "Control")$CLR
clr_eQTL <- filter(north_clr_data_long, Method == "eQTL")$CLR
wilcox.test(clr_control, clr_lfmm) # p-value < 2.2e-16
wilcox.test(clr_control, clr_pcad) # p-value = 6.97e-07
wilcox.test(clr_pcad, clr_lfmm) # p-value < 2.2e-16
wilcox.test(clr_control, clr_eQTL)


# Plot CLR distributions for Coastal populations
coast_clr_LFMM <- coast_lfmm$LR
coast_clr_PCAdapt <- coast_pcadapt$LR
coast_clr_eQTL <- coast_eQTL$LR
coast_clr_eGene <- coast_eGenes$LR
coast_clr_control <- coast_control$LR

coast_clr_data <- cbind(coast_clr_LFMM, coast_clr_PCAdapt, coast_clr_eGene, coast_clr_eQTL, coast_clr_control)

coast_clr_data <- as.data.frame(coast_clr_data) %>%
  na.omit()

names(coast_clr_data) <- c("LFMM", "PCAdapt", "eQTL", "eGene", "Control")

coast_clr_data_long <- coast_clr_data %>%
  gather(Method, CLR) %>%
  mutate(CLR = log10(CLR + 1))

# Change CLR column to numeric
coast_clr_data_long[, 2] <- data.frame(sapply(coast_clr_data_long[, 2], as.numeric))

coast_clr_data_long <- coast_clr_data_long[is.finite(coast_clr_data_long$CLR), ] %>%
  na.omit()

median(filter(coast_clr_data_long, Method == "Control")$CLR) # 0.02097955

coast_clr_lfmm <- filter(coast_clr_data_long, Method == "LFMM")$CLR
coast_clr_pcad <- filter(coast_clr_data_long, Method == "PCAdapt")$CLR
coast_clr_control <- filter(coast_clr_data_long, Method == "Control")$CLR
wilcox.test(coast_clr_control, coast_clr_lfmm) # p-value = 0.05523
wilcox.test(coast_clr_control, coast_clr_pcad) # p-value < 2.2e-16
wilcox.test(coast_clr_pcad, coast_clr_lfmm) # p-value < 2.2e-16

p2 <- coast_clr_data_long %>%
  ggplot(aes(x=Method, y=CLR, fill=Method)) +
  geom_violin(adjust = 13, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.25, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.054, color="black", alpha = 0.8) +
  scale_y_continuous(trans = 'log10') +
  xlab(" ") +
  ylab("CLR") +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Coast_CLR.tiff", plot = p2, device = "tiff")

clr_plot <- plot_grid(p, p2)
ggsave("PLOS_genetics_figures/North_Coast_CLR.tiff", plot = clr_plot, device = "tiff", width = 10, height = 5)


# Plot CLR for selection outliers and shared outliers between LFMM and PCAdapt (coastal populations)
coast_clr_LFMM <- coast_lfmm$LR
coast_clr_PCAdapt <- coast_pcadapt$LR
coast_clr_shared <- coast_shared$LR

coast_clr_data <- cbind(coast_clr_LFMM, coast_clr_PCAdapt, coast_clr_shared)

coast_clr_data <- as.data.frame(coast_clr_data) %>%
  na.omit()

names(coast_clr_data) <- c("LFMM", "PCAdapt", "Shared")

coast_clr_data_long <- coast_clr_data %>%
  gather(Method, CLR) %>%
  mutate(CLR = log10(CLR + 1))

# Change CLR column to numeric
coast_clr_data_long[, 2] <- data.frame(sapply(coast_clr_data_long[, 2], as.numeric))

coast_clr_data_long <- coast_clr_data_long[is.finite(coast_clr_data_long$CLR), ] %>%
  na.omit()

median(filter(coast_clr_data_long, Method == "Shared")$CLR) # 0.01184523

coast_clr_lfmm <- filter(coast_clr_data_long, Method == "LFMM")$CLR
coast_clr_pcad <- filter(coast_clr_data_long, Method == "PCAdapt")$CLR
coast_clr_control <- filter(coast_clr_data_long, Method == "Shared")$CLR
wilcox.test(coast_clr_control, coast_clr_lfmm) # p-value = 0.05523
wilcox.test(coast_clr_shared, coast_clr_pcad) # p-value < 2.2e-16
wilcox.test(coast_clr_pcad, coast_clr_lfmm) # p-value < 2.2e-16

# Define color palette for groups
my_colors3 <- c("Shared" = "darkgrey", "LFMM" = "aquamarine4", "PCAdapt" = "brown")

# Define the plot order
coast_clr_data_long$Method <- factor(coast_clr_data_long$Method , levels = c("Shared", "LFMM", "PCAdapt"))

q1 <- coast_clr_data_long %>%
  ggplot(aes(x=Method, y=CLR, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.25, label = "p.signif",
                     ref.group = "Shared") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.113, color="black", alpha = 0.8) +
  scale_y_continuous(trans = 'log10') +
  xlab(" ") +
  ylab("CLR") +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors3) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Coast_shared_CLR.tiff", plot = q1, device = "tiff")


# Plot CLR for selection outliers and shared outliers between LFMM and PCAdapt (Northern populations)
north_clr_LFMM <- north_lfmm$LR
north_clr_PCAdapt <- north_pcadapt$LR
north_clr_shared <- north_shared$LR

north_clr_data <- cbind(north_clr_LFMM, north_clr_PCAdapt, north_clr_shared)

north_clr_data <- as.data.frame(north_clr_data) %>%
  na.omit()

names(north_clr_data) <- c("LFMM", "PCAdapt", "Shared")

north_clr_data_long <- north_clr_data %>%
  gather(Method, CLR) %>%
  mutate(CLR = log10(CLR + 1))

# Change CLR column to numeric
north_clr_data_long[, 2] <- data.frame(sapply(north_clr_data_long[, 2], as.numeric))

north_clr_data_long <- north_clr_data_long[is.finite(north_clr_data_long$CLR), ] %>%
  na.omit()

median(filter(north_clr_data_long, Method == "Shared")$CLR) # 0

north_clr_lfmm <- filter(north_clr_data_long, Method == "LFMM")$CLR
north_clr_pcad <- filter(north_clr_data_long, Method == "PCAdapt")$CLR
north_clr_control <- filter(north_clr_data_long, Method == "Shared")$CLR
wilcox.test(north_clr_shared, north_clr_lfmm) # p-value = 0.05523
wilcox.test(north_clr_shared, north_clr_pcad) # p-value < 2.2e-16
wilcox.test(north_clr_pcad, north_clr_lfmm) # p-value < 2.2e-16

# Define color palette for groups
my_colors3 <- c("Shared" = "darkgrey", "LFMM" = "aquamarine4", "PCAdapt" = "brown")

# Define the plot order
north_clr_data_long$Method <- factor(north_clr_data_long$Method , levels = c("Shared", "LFMM", "PCAdapt"))

q2 <- north_clr_data_long %>%
  ggplot(aes(x=Method, y=CLR, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.25, label = "p.signif",
                     ref.group = "Shared") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.113, color="black", alpha = 0.8) +
  scale_y_continuous(trans = 'log10') +
  xlab(" ") +
  ylab("CLR") +
  ggtitle("North") +
  scale_fill_manual(values = my_colors3) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_shared_CLR.tiff", plot = q2, device = "tiff")

shared_clr_plot <- plot_grid(q2, q1)
ggsave("PLOS_genetics_figures/shared_clr_plot.tiff", plot = shared_clr_plot, device = "tiff", width = 10, height = 5)

#------------------------------------------------------------------------

# Plot Fay and Wu's H for North population 
north_fayh_LFMM <-north_lfmm_theta$fayH
north_fayh_PCAdapt <- north_PCAdapt_theta$fayH
north_fayh_eQTL <- north_eQTL_theta$fayH
north_fayh_eGene <- north_eGenes_theta$fayH
north_fayh_Control <- north_control_theta$fayH
north_fayh_data <- cbind(north_fayh_LFMM, north_fayh_PCAdapt, north_fayh_eGene, north_fayh_eQTL, north_fayh_Control)

north_fayh_data <- as.data.frame(north_fayh_data) %>%
  na.omit()

names(north_fayh_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

north_fayh_data_long <- north_fayh_data %>%
  gather(Method, FayH)

# Change CLR column to numeric
north_fayh_data_long[, 2] <- data.frame(sapply(north_fayh_data_long[, 2], as.numeric)) 

north_fayh_data_long <- north_fayh_data_long[is.finite(north_fayh_data_long$FayH), ] %>%
  na.omit() 

median(filter(north_fayh_data_long, Method == "Control")$FayH) # 0.030783

north_fay_lfmm <- filter(north_fayh_data_long, Method == "LFMM")$FayH
north_fay_pcad <- filter(north_fayh_data_long, Method == "PCAdapt")$FayH
north_fay_control <- filter(north_fayh_data_long, Method == "Control")$FayH
wilcox.test(north_fay_control, north_fay_lfmm) # p-value = 0.001
wilcox.test(north_fay_control, north_fay_pcad) # p-value = 6.623e-06
wilcox.test(north_fay_pcad, north_fay_lfmm) # p-value = 0.01233

p3 <- north_fayh_data_long %>%
  ggplot(aes(x=Method, y=FayH, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.8, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.030783, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(expression(paste("Fay & Wu's ", italic("H")))) +
  ggtitle("North") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_FayH.tiff", plot = p3, device = "tiff")


# Plot FayH for coast population
coast_fayh_LFMM <- coast_lfmm_theta$fayH
coast_fayh_PCAdapt <- coast_PCAdapt_theta$fayH
coast_fayh_eQTL <- coast_eQTL_theta$fayH
coast_fayh_eGene <- coast_eGenes_theta$fayH
coast_fayh_Control <- coast_control_theta$fayH

coast_fayh_data <- cbind(coast_fayh_LFMM, coast_fayh_PCAdapt, coast_fayh_eGene, coast_fayh_eQTL, coast_fayh_Control)

coast_fayh_data <- as.data.frame(coast_fayh_data) %>%
  na.omit()

names(coast_fayh_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

coast_fayh_data_long <- coast_fayh_data %>%
  gather(Method, FayH) 

coast_fayh_data_long <- coast_fayh_data_long[is.finite(coast_fayh_data_long$FayH), ]

median(filter(coast_fayh_data_long, Method == "Control")$FayH) # 0.216144

coast_fayh_lfmm <- filter(coast_fayh_data_long, Method == "LFMM")$FayH
coast_fayh_pcad <- filter(coast_fayh_data_long, Method == "PCAdapt")$FayH
coast_fayh_control <- filter(coast_fayh_data_long, Method == "Control")$FayH
wilcox.test(coast_fayh_control, coast_fayh_lfmm) # p-value = 4.998e-05
wilcox.test(coast_fayh_control, coast_fayh_pcad) # p-value = 5.977e-05
wilcox.test(coast_fayh_pcad, coast_fayh_lfmm) # p-value = 0.6726

p4 <- coast_fayh_data_long %>%
  ggplot(aes(x=Method, y=FayH, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2.4, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.216144, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(expression(paste("Fay & Wu's ", italic("H")))) +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Coast_FayH.tiff", plot = p4, device = "tiff")

FayH_plot <- plot_grid(p3, p4)
ggsave("PLOS_genetics_figures/North_Coast_FayH.tiff", plot = FayH_plot, device = "tiff", width = 10, height = 5)


# Plot Shared Fay and Wu's H for North population 
north_fayh_LFMM <-north_lfmm_theta$fayH
north_fayh_PCAdapt <- north_PCAdapt_theta$fayH
north_fayh_shared <- north_shared_theta$fayH

north_fayh_data <- cbind(north_fayh_LFMM, north_fayh_PCAdapt, north_fayh_shared)

north_fayh_data <- as.data.frame(north_fayh_data) %>%
  na.omit()

names(north_fayh_data) <- c("LFMM", "PCAdapt", "Shared")

north_fayh_data_long <- north_fayh_data %>%
  gather(Method, FayH)

# Change CLR column to numeric
north_fayh_data_long[, 2] <- data.frame(sapply(north_fayh_data_long[, 2], as.numeric)) 

north_fayh_data_long <- north_fayh_data_long[is.finite(north_fayh_data_long$FayH), ] %>%
  na.omit() 

median(filter(north_fayh_data_long, Method == "Shared")$FayH) # 0.254183

north_fay_lfmm <- filter(north_fayh_data_long, Method == "LFMM")$FayH
north_fay_pcad <- filter(north_fayh_data_long, Method == "PCAdapt")$FayH
north_fay_shared <- filter(north_fayh_data_long, Method == "Shared")$FayH
wilcox.test(north_fay_shared, north_fay_lfmm) # p-value = 0.0001314
wilcox.test(north_fay_control, north_fay_pcad) # p-value = 6.623e-06
wilcox.test(north_fay_pcad, north_fay_lfmm) # p-value = 0.01233

# Define the plot order
north_fayh_data_long$Method <- factor(north_fayh_data_long$Method , levels = c("Shared", "LFMM", "PCAdapt"))

q3 <- north_fayh_data_long %>%
  ggplot(aes(x=Method, y=FayH, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2, label = "p.signif",
                     ref.group = "Shared") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.262442, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(expression(paste("Fay & Wu's ", italic("H")))) +
  ggtitle("North") +
  scale_fill_manual(values = my_colors3) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_Shared_FayH.tiff", plot = q3, device = "tiff")


# Plot Shared Fay and Wu's H for Coast population 
coast_fayh_LFMM <- coast_lfmm_theta$fayH
coast_fayh_PCAdapt <- coast_PCAdapt_theta$fayH
coast_fayh_shared <- coast_shared_theta$fayH

coast_fayh_data <- cbind(coast_fayh_LFMM, coast_fayh_PCAdapt, coast_fayh_shared)

coast_fayh_data <- as.data.frame(coast_fayh_data) %>%
  na.omit()

names(coast_fayh_data) <- c("LFMM", "PCAdapt", "Shared")

coast_fayh_data_long <- coast_fayh_data %>%
  gather(Method, FayH)

# Change CLR column to numeric
coast_fayh_data_long[, 2] <- data.frame(sapply(coast_fayh_data_long[, 2], as.numeric)) 

coast_fayh_data_long <- coast_fayh_data_long[is.finite(coast_fayh_data_long$FayH), ] %>%
  na.omit() 

median(filter(coast_fayh_data_long, Method == "Shared")$FayH) # 0.362442

coast_fay_lfmm <- filter(coast_fayh_data_long, Method == "LFMM")$FayH
coast_fay_pcad <- filter(coast_fayh_data_long, Method == "PCAdapt")$FayH
coast_fay_shared <- filter(coast_fayh_data_long, Method == "Shared")$FayH
wilcox.test(coast_fay_shared, coast_fay_lfmm) # p-value = 0.01083
wilcox.test(coast_fay_shared, coast_fay_pcad) # p-value = 0.0005292
wilcox.test(coast_fay_pcad, coast_fay_lfmm) # p-value = 0.1098

# Define the plot order
coast_fayh_data_long$Method <- factor(coast_fayh_data_long$Method , levels = c("Shared", "LFMM", "PCAdapt"))

q4 <- coast_fayh_data_long %>%
  ggplot(aes(x=Method, y=FayH, fill=Method)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 2, label = "p.signif",
                     ref.group = "Shared") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.262442, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(expression(paste("Fay & Wu's ", italic("H")))) +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors3) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/coast_Shared_FayH.tiff", plot = q4, device = "tiff")

Shared_FayH_plot <- plot_grid(q3, q4)
ggsave("PLOS_genetics_figures/Shared_FayH_plot.tiff", plot = Shared_FayH_plot, device = "tiff", width = 10, height = 5)


#---------------------------------------------------------------------
# Plot Watterson thetas for North population
north_wat_theta_LFMM <- north_lfmm_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

north_wat_theta_LFMM <- north_wat_theta_LFMM$Watterson_theta

north_wat_theta_PCAdapt <- north_PCAdapt_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

north_wat_theta_PCAdapt <- north_wat_theta_PCAdapt$Watterson_theta

north_wat_theta_eGene <- north_eGenes_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

north_wat_theta_eGene <- north_wat_theta_eGene$Watterson_theta

north_wat_theta_eQTL <- north_eQTL_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

north_wat_theta_eQTL <- north_wat_theta_eQTL$Watterson_theta

north_wat_theta_Control <- north_control_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

north_wat_theta_Control <- north_wat_theta_Control$Watterson_theta

north_wat_theta_data <- cbind(north_wat_theta_LFMM, north_wat_theta_PCAdapt, north_wat_theta_eGene, north_wat_theta_eQTL, north_wat_theta_Control)

north_wat_theta_data <- as.data.frame(north_wat_theta_data)
names(north_wat_theta_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

north_wat_theta_data_long <- north_wat_theta_data %>%
  gather(Method, Watterson_theta) %>%
  na.omit()

north_wat_theta_data_long <- north_wat_theta_data_long[is.finite(north_wat_theta_data_long$Watterson_theta), ]

median(filter(north_wat_theta_data_long, Method == "Control")$Watterson_theta) # 0.2200703

north_wat_theta_lfmm <- filter(north_wat_theta_data_long, Method == "LFMM")$Watterson_theta
north_wat_theta_pcad <- filter(north_wat_theta_data_long, Method == "PCAdapt")$Watterson_theta
north_wat_theta_control <- filter(north_wat_theta_data_long, Method == "Control")$Watterson_theta
wilcox.test(north_wat_theta_control, north_wat_theta_lfmm) # p-value = 0.06571
wilcox.test(north_wat_theta_control, north_wat_theta_pcad) # p-value = 0.003819
wilcox.test(north_wat_theta_pcad, north_wat_theta_lfmm) # p-value = p-value < 2.2e-16

p5 <- north_wat_theta_data_long %>%
  ggplot(aes(x=Method, y=Watterson_theta, fill=(Method))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 0.4, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  geom_hline(yintercept = 0.2200703, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(bquote("θ" ~ italic("w"))) +
  ggtitle("North") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_watterson_theta.tiff", plot = p5, device = "tiff")


#------------------------------------------------------------------------------
# Plot watterson thetas for Coast population
coast_wat_theta_LFMM <- coast_lfmm_theta %>%
  dplyr::select(., tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

coast_wat_theta_LFMM <- coast_wat_theta_LFMM$Watterson_theta

coast_wat_theta_PCAdapt <- coast_PCAdapt_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

coast_wat_theta_PCAdapt <- coast_wat_theta_PCAdapt$Watterson_theta

coast_wat_theta_eGene <- coast_eGenes_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

coast_wat_theta_eGene <- coast_wat_theta_eGene$Watterson_theta

coast_wat_theta_eQTL <- coast_eQTL_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

coast_wat_theta_eQTL <- coast_wat_theta_eQTL$Watterson_theta

coast_wat_theta_Control <- coast_control_theta %>%
  dplyr::select(, tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites)

coast_wat_theta_Control <- coast_wat_theta_Control$Watterson_theta

coast_wat_theta_data <- cbind(coast_wat_theta_LFMM, coast_wat_theta_PCAdapt, coast_wat_theta_eGene, coast_wat_theta_eQTL, coast_wat_theta_Control)

coast_wat_theta_data <- as.data.frame(coast_wat_theta_data)

names(coast_wat_theta_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

coast_wat_theta_data_long <- coast_wat_theta_data %>%
  gather(Method, Watterson_theta) 

coast_wat_theta_data_long <- coast_wat_theta_data_long[is.finite(coast_wat_theta_data_long$Watterson_theta), ]

median(filter(coast_wat_theta_data_long, Method == "Control")$Watterson_theta) # 0.2881245

coast_wat_theta_lfmm <- filter(coast_wat_theta_data_long, Method == "LFMM")$Watterson_theta
coast_wat_theta_pcad <- filter(coast_wat_theta_data_long, Method == "PCAdapt")$Watterson_theta
coast_wat_theta_control <- filter(coast_wat_theta_data_long, Method == "Control")$Watterson_theta
wilcox.test(coast_wat_theta_control, coast_wat_theta_lfmm) # p-value = 1.57e-14
wilcox.test(coast_wat_theta_control, coast_wat_theta_pcad) # p-value = 0.0008722
wilcox.test(coast_wat_theta_pcad, coast_wat_theta_lfmm) # p-value < 2.2e-16

p6 <- coast_wat_theta_data_long %>%
  ggplot(aes(x=Method, y=Watterson_theta, fill=(Method))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 0.4, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference            
  geom_hline(yintercept = 0.2881245, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(bquote("θ" ~ italic("w"))) +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Coast_watterson_theta.tiff", plot = p6, device = "tiff")

North_Coast_watterson_theta_plot <- plot_grid(p5, p6)
ggsave("PLOS_genetics_figures/North_Coast_watterson_theta_plot.tiff", plot = North_Coast_watterson_theta_plot, device = "tiff", width = 10, height = 5)

#------------------------------------------------------------------------------------
# Plot thetas pie for North population
north_theta_pie_LFMM <- north_lfmm_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

north_theta_pie_LFMM <- north_theta_pie_LFMM$theta_pie

north_theta_pie_PCAdapt <- north_PCAdapt_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

north_theta_pie_PCAdapt <- north_theta_pie_PCAdapt$theta_pie

north_theta_pie_eGene <- north_eGenes_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

north_theta_pie_eGene <- north_theta_pie_eGene$theta_pie

north_theta_pie_eQTL <- north_eQTL_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

north_theta_pie_eQTL <- north_theta_pie_eQTL$theta_pie

north_theta_pie_Control <- north_control_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

north_theta_pie_Control <- north_theta_pie_Control$theta_pie

north_theta_pie_data <- cbind(north_theta_pie_LFMM, north_theta_pie_PCAdapt, north_theta_pie_eGene, north_theta_pie_eQTL, north_theta_pie_Control)

north_theta_pie_data <- as.data.frame(north_theta_pie_data)
names(north_theta_pie_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

north_theta_pie_data_long <- north_theta_pie_data %>%
  gather(Method, theta_pie) 

north_theta_pie_data_long <- north_theta_pie_data_long[is.finite(north_theta_pie_data_long$theta_pie), ]

median(filter(north_theta_pie_data_long, Method == "Control")$theta_pie) # 0.1834275

# Wilcoxon test
north_theta_pie_lfmm <- filter(north_theta_pie_data_long, Method == "LFMM")$theta_pie
north_theta_pie_pcad <- filter(north_theta_pie_data_long, Method == "PCAdapt")$theta_pie
north_theta_pie_control <- filter(north_theta_pie_data_long, Method == "Control")$theta_pie
wilcox.test(north_theta_pie_control, north_theta_pie_lfmm) # p-value < 2.2e-16
wilcox.test(north_theta_pie_control, north_theta_pie_pcad) # p-value = 0.3752
wilcox.test(north_theta_pie_pcad, north_theta_pie_lfmm) # p-value < 2.2e-16


p7 <- north_theta_pie_data_long %>%
  ggplot(aes(x=Method, y=theta_pie, fill=(Method))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 0.7, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference    
  geom_hline(yintercept = 0.1834275, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(bquote("θ" ~ italic("π"))) +
  ggtitle("North") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/North_theta_pie.tiff", plot = p7, device = "tiff")

#---------------------------------------------------------------------------------

# Plot thetas pie for coast population
coast_theta_pie_LFMM <- coast_lfmm_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

coast_theta_pie_LFMM <- coast_theta_pie_LFMM$theta_pie

coast_theta_pie_PCAdapt <- coast_PCAdapt_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

coast_theta_pie_PCAdapt <- coast_theta_pie_PCAdapt$theta_pie

coast_theta_pie_eGene <- coast_eGenes_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

coast_theta_pie_eGene <- coast_theta_pie_eGene$theta_pie

coast_theta_pie_eQTL <- coast_eQTL_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

coast_theta_pie_eQTL <- coast_theta_pie_eQTL$theta_pie

coast_theta_pie_Control <- coast_control_theta %>%
  dplyr::select(, tP, numSites) %>%
  mutate(theta_pie = tP/numSites)

coast_theta_pie_Control <- coast_theta_pie_Control$theta_pie

coast_theta_pie_data <- cbind(coast_theta_pie_LFMM, coast_theta_pie_PCAdapt, coast_theta_pie_eGene, coast_theta_pie_eQTL, coast_theta_pie_Control)

coast_theta_pie_data <- as.data.frame(coast_theta_pie_data)

names(coast_theta_pie_data) <- c("LFMM", "PCAdapt", "eGene", "eQTL", "Control")

coast_theta_pie_data_long <- coast_theta_pie_data %>%
  gather(Method, theta_pie) 

coast_theta_pie_data_long <- coast_theta_pie_data_long[is.finite(coast_theta_pie_data_long$theta_pie), ]

median(filter(coast_theta_pie_data_long, Method == "Control")$theta_pie) # 0.2172999

p8 <- coast_theta_pie_data_long %>%
  ggplot(aes(x=Method, y=theta_pie, fill=(Method))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 0.7, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference    
  geom_hline(yintercept = 0.2172999, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab(bquote("θ" ~ italic("π"))) +
  ggtitle("Coast") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Coast_theta_pie.tiff", plot = p8, device = "tiff")

# Wilcoxon test
coast_theta_pie_lfmm <- filter(coast_theta_pie_data_long, Method == "LFMM")$theta_pie
coast_theta_pie_pcad <- filter(coast_theta_pie_data_long, Method == "PCAdapt")$theta_pie
coast_theta_pie_control <- filter(coast_theta_pie_data_long, Method == "Control")$theta_pie
wilcox.test(coast_theta_pie_control, coast_theta_pie_lfmm) # p-value = 0.2015
wilcox.test(coast_theta_pie_control, coast_theta_pie_pcad) # p-value = 0.004505
wilcox.test(coast_theta_pie_pcad, coast_theta_pie_lfmm) # p-value < 2.2e-16


North_Coast_theta_pie_plot <- plot_grid(p7, p8)
ggsave("PLOS_genetics_figures/North_Coast_theta_pie_plot.tiff", plot = North_Coast_theta_pie_plot, device = "tiff", width = 10, height = 5)

#------------------------------------------------------------------------------

# Plot FSTs for LFMM, PCAdapt, and Control
LFMM <- LFMM_FST$weighted_FST
PCAdapt <- PCAdapt_FST$weighted_FST
eQTL <- eQTL_FST$weighted_FST
eGene <- eGene_FST$weighted_FST
Control <- Control_FST$weighted_FST

data9 <- cbind(LFMM, PCAdapt, eGene, eQTL, Control)

data9 <- as.data.frame(data9)

data9_long <- data9 %>%
  gather(Method, FST)

data9_long <- data9_long[is.finite(data9_long$FST), ]

median(filter(data9_long, Method == "Control")$FST) # 0.164796

p9 <- data9_long %>%
  ggplot(aes(x=Method, y=FST, fill=(Method))) +
  geom_violin(adjust = 10, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  # stat_summary(fun=mean, geom="point", shape=20, size=4) +
  stat_compare_means(method = "wilcox.test", label.y = 0.8, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference    
  geom_hline(yintercept = 0.164796, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab("FST") +
  ggtitle(" ") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/FST.tiff", plot = p9, device = "tiff", width = 7, height = 4)


FST_lfmm <- filter(data9_long, Method == "LFMM")$FST
FST_pcad <- filter(data9_long, Method == "PCAdapt")$FST
FST_control <- filter(data9_long, Method == "Control")$FST
wilcox.test(FST_control, FST_lfmm) # p-value < 2.2e-16
wilcox.test(FST_control, FST_pcad) # p-value = 3.941e-06
wilcox.test(FST_pcad, FST_lfmm) # p-value < 2.2e-16

#-------------------------------------------------------------------------------
# Plot Dxy 
LFMM <- LFMM_dxy$dxy_stat
PCAdapt <- PCAdapt_dxy$dxy_stat
eQTL <- eQTL_dxy$dxy_stat
eGene <- eGene_dxy$dxy_stat
Control <- Control_dxy$dxy_stat

data10 <- cbind(LFMM, PCAdapt, eGene, eQTL, Control)

data10 <- as.data.frame(data10)
#data10 <- log(data10)

data10_long <- data10 %>%
  gather(Method, DXY) 

data10_long <- data10_long[is.finite(data10_long$DXY), ]

median(filter(data10_long, Method == "Control")$DXY) # 0.006734933

DXY_lfmm <- filter(data10_long, Method == "LFMM")$DXY
DXY_pcad <- filter(data10_long, Method == "PCAdapt")$DXY
DXY_control <- filter(data10_long, Method == "Control")$DXY
wilcox.test(DXY_control, DXY_lfmm) # p-value < 2.2e-16
wilcox.test(DXY_control, DXY_pcad) # p-value < 2.2e-16
wilcox.test(DXY_pcad, DXY_lfmm) # p-value = 0.8483

p10 <- data10_long %>%
  ggplot(aes(x=Method, y=DXY, fill=(Method))) +
  geom_violin(adjust = 10, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label.y = 0.3, label = "p.signif",
                     ref.group = "Control") +      # Add global p-value
  # Pairwise comparison against reference                # Pairwise comparison against reference
  scale_y_continuous(trans = 'log10') +
  geom_hline(yintercept = 0.006734933, color="black", alpha = 0.8) +
  xlab(" ") +
  ylab("Dxy") +
  ggtitle(" ") +
  scale_fill_manual(values = my_colors) +  # Assign colors to groups
  guides(fill = "none") +
  theme_bw() 

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/Dxy.tiff", plot = p10, device = "tiff", width = 7, height = 4)

FST_Dxy_plot <- plot_grid(p9, p10)
ggsave("PLOS_genetics_figures/FST_Dxy_plot.tiff", plot = FST_Dxy_plot, device = "tiff", width = 10, height = 5)

#---------------------------------------------------------------
lfmm_genes <- read.table("../argo_plots/lfmm_genes.txt") %>%
  dplyr::select(V3) %>%
  separate(V3, into = c("name", "gene")) %>%
  filter(!gene %in% c("", "prime")) %>%
  dplyr::select(gene) %>%
  distinct() # 560 genes housing lfmm outliers

pcadapt_genes <- read.table("../argo_plots/pcadapt_genes.txt") %>%
  dplyr::select(V3) %>%
  separate(V3, into = c("name", "gene")) %>%
  filter(!gene %in% c("", "prime")) %>%
  dplyr::select(gene) %>%
  distinct() # 1122 genes housing PCAdapt outliers

eQTL_genes <- read.table("../network_analysis/eQTL_genes2.txt") %>%
  #dplyr::select(V3) %>%
  #separate(V3, into = c("name", "gene")) %>%
  #filter(!gene %in% c("", "prime")) %>%
  #dplyr::select(gene) %>%
  filter_all(all_vars(. != ".")) %>%
  distinct() #1824 eQTL genes

eGene_genes <- read.table("../eQTL_analysis/eGenes.txt")  %>%
  #dplyr::select(V3) %>%
  #separate(V3, into = c("name", "gene")) %>%
  #filter(!gene %in% c("", "prime")) %>%
  #dplyr::select(gene) %>%
  filter_all(all_vars(. != ".")) %>%
  distinct() # 3129 genes are eGenes

All_genes <- read.table("../network_analysis/all_genes_eQTL.txt") %>%
  #dplyr::select(V3) %>%
  #separate(V3, into = c("name", "gene")) %>%
  #filter(!gene %in% c("", "prime")) %>%
  #dplyr::select(gene) %>%
  filter_all(all_vars(. != ".")) %>%
  distinct() # Total number of genes investigated are 42638


LFMM <- intersect(eQTL_genes$V1, lfmm_genes$gene) #114/560; 20.4% of LFMM genes are eQTLs
PCAdapt <- intersect(eQTL_genes$V1, pcadapt_genes$gene) #208/1122 18.5% of PCAdapt genes are eQTLs
All <- intersect(eQTL_genes$V1, All_genes$Gene_ID) #5734/42638 13.4% of all genes are eQTls. 


eGene_LFMM <- intersect(eGene_genes$V1, lfmm_genes$gene) # 11/560; 2% of lfmm genes are eGenes
eGene_pcadapt <- intersect(eGene_genes$V1, pcadapt_genes$gene) # 24/1122 2.1% of PCAdapt genes are eGenes
eGene_All <- intersect(eGene_genes$V1 , All_genes$Gene_ID) # 3128/42638 7.3%

eGene_eQTL <- intersect(eGene_genes$V1, eQTL_genes$V1) #26/42638 (0.06%) are both eQTLs and eGenes. 


# Get the number of times eQTLs are found in LFMM, PCAdapt, and across all genes
LFMM_count <- tibble(letters = lfmm_genes$gene, count = unlist(map(lfmm_genes$gene, function(x) sum(eQTL_genes$V1 %in% x))))
colnames(LFMM_count)[2] <- "LFMM"

PCAdapt_count <- tibble(letters = pcadapt_genes$gene, count = unlist(map(pcadapt_genes$gene, function(x) sum(eQTL_genes$V1 %in% x))))
colnames(PCAdapt_count)[2] <- "PCAdapt"

All_count <- tibble(letters = All_genes$Gene_ID, count = unlist(map(All_genes$Gene_ID, function(x) sum(eQTL_genes$V1 %in% x))))
colnames(All_count)[2] <- "All"

# Get the number of times eGenes are found in LFMM, PCAdapt, and across all genes
LFMM2_count <- tibble(letters = lfmm_genes$gene, count = unlist(map(lfmm_genes$gene, function(x) sum(eGene_genes$V1 %in% x))))
colnames(LFMM_count)[2] <- "LFMM"

PCAdapt2_count <- tibble(letters = pcadapt_genes$gene, count = unlist(map(pcadapt_genes$gene, function(x) sum(eGene_genes$V1 %in% x))))
colnames(PCAdapt_count)[2] <- "PCAdapt"

All2_count <- tibble(letters = All_genes$Gene_ID, count = unlist(map(All_genes$Gene_ID, function(x) sum(eGene_genes$V1 %in% x))))
colnames(All_count)[2] <- "All"

# Add all eQTL occurences in one dataframe
percent_eQTL <- cbind(LFMM_count$LFMM, PCAdapt_count$PCAdapt, All_count$All)

percent_eQTL <- data.frame(percent_eQTL)
names(percent_eQTL) <- c("LFMM", "PCAdapt", "All")

percent_eQTL <- percent_eQTL %>%
  gather(eQTL, count) 

# Add all eGene occurences in one dataframe
percent_eGene <- cbind(LFMM2_count$count, PCAdapt2_count$count, All2_count$count)

percent_eGene <- data.frame(percent_eGene)
names(percent_eGene) <- c("LFMM", "PCAdapt", "All")

percent_eGene <- percent_eGene %>%
  gather(eGene, count) 

# Assign new colors
my_colors2 <- c("All" = "darkgrey", "LFMM" = "brown", "PCAdapt" = "aquamarine4")

# Get eQTL proportions
percent_eQTL2 <- percent_eQTL %>%
  group_by(eQTL) %>%
  summarise(total = n(),
            successes = sum(count == 1)) %>%
  mutate(proportion = successes / total)


# Get eGene proportions
percent_eGene2 <- percent_eGene %>%
  group_by(eGene) %>%
  summarise(total = n(),
            successes = sum(count == 1)) %>%
  mutate(proportion = successes / total)

# Calculate proportions and Fisher's exact test
eQTL_summary <- percent_eQTL2 %>%
  group_by(eQTL) %>%
  summarise(proportion = mean(proportion),
            total = mean(total)) %>%
  ungroup()

eQTL_summary <- eQTL_summary %>%
  mutate(count = proportion * total,
         non_count = total - count)

# Fisher's exact test calculations
pcadapt_total <- eQTL_summary$total[eQTL_summary$eQTL == "PCAdapt"]
all_total <- eQTL_summary$total[eQTL_summary$eQTL == "All"]
pcadapt_count <- eQTL_summary$count[eQTL_summary$eQTL == "PCAdapt"]
all_count <- eQTL_summary$count[eQTL_summary$eQTL == "All"]

# Construct contingency table
contingency_table <- matrix(c(
  pcadapt_count, pcadapt_total - pcadapt_count,
  all_count, all_total - all_count
), nrow = 2)

fisher_test <- fisher.test(contingency_table)
p_value <- fisher_test$p.value

# Bootstrap and plot
p11 <- data.frame(boot = 1:1000) %>%
  group_by(boot) %>% 
  do(sample_n(percent_eQTL2, nrow(percent_eQTL2), replace = TRUE)) %>%
  group_by(boot, eQTL) %>%
  mutate(ci_low = proportion - qnorm(0.975) * sqrt(proportion * (1 - proportion) / total),
         ci_high = proportion + qnorm(0.975) * sqrt(proportion * (1 - proportion) / total)) %>%
  group_by(eQTL) %>% 
  summarise(proportion = mean(proportion),
            ci_low = mean(ci_low),
            ci_high = mean(ci_high)) %>%
  ggplot(aes(x = eQTL, y = proportion, fill = eQTL)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, position = position_dodge(0.9)) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("eQTL") +
  ylab("Percentage") +
  ggtitle(" ") +
  scale_fill_manual(values = my_colors2) +  # Adjust to your color preferences
  guides(fill = "none") +
  theme_bw() +
  annotate("text", x = 1.5, y = max(percent_eQTL2$proportion) + 0.02, 
           label = paste("Fisher's p-value:", format("< 2.2e-16", digits = 3)), size = 4, hjust = 0.5)


# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/eQTL_proportions.tiff", plot = p11, device = "tiff", width = 10, height = 10)


# Get eGene proportions
percent_eGene2 <- percent_eGene %>%
  group_by(eGene) %>%
  summarise(total = n(),
            successes = sum(count == 1)) %>%
  mutate(proportion = successes / total)

# Calculate proportions and Fisher's exact test
eGene_summary <- percent_eGene2 %>%
  group_by(eGene) %>%
  summarise(proportion = mean(proportion),
            total = mean(total)) %>%
  ungroup()

eGene_summary <- eGene_summary %>%
  mutate(count = proportion * total,
         non_count = total - count)

# Fisher's exact test calculations
pcadapt_total <- eGene_summary$total[eGene_summary$eGene == "PCAdapt"]
all_total <- eGene_summary$total[eGene_summary$eGene == "All"]
pcadapt_count <- eGene_summary$count[eGene_summary$eGene == "PCAdapt"]
all_count <- eGene_summary$count[eGene_summary$eGene == "All"]

# Construct contingency table
contingency_table <- matrix(c(
  pcadapt_count, pcadapt_total - pcadapt_count,
  all_count, all_total - all_count
), nrow = 2)

fisher_test <- fisher.test(contingency_table)
p_value <- fisher_test$p.value

## No significant difference between LFMM and All eGenes proportions; p = 0.632

# Bootstrap and plot
p12 <- data.frame(boot = 1:1000) %>%
  group_by(boot) %>% 
  do(sample_n(percent_eGene2, nrow(percent_eGene2), replace = TRUE)) %>%
  group_by(boot, eGene) %>%
  mutate(ci_low = proportion - qnorm(0.975) * sqrt(proportion * (1 - proportion) / total),
         ci_high = proportion + qnorm(0.975) * sqrt(proportion * (1 - proportion) / total)) %>%
  group_by(eGene) %>% 
  summarise(proportion = mean(proportion),
            ci_low = mean(ci_low),
            ci_high = mean(ci_high)) %>%
  ggplot(aes(x = eGene, y = proportion, fill = eGene)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, position = position_dodge(0.9)) +
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("eGene") +
  ylab("Percentage") +
  ggtitle(" ") +
  scale_fill_manual(values = my_colors2) +  # Adjust to your color preferences
  guides(fill = "none") +
  theme_bw() +
  annotate("text", x = 1.5, y = max(percent_eGene2$proportion) + 0.09, 
           label = paste("Fisher's p-value:", format("< 2.2e-16", digits = 3)), size = 4, hjust = 0.5)

# Save the ggplot as a .tiff file
ggsave("PLOS_genetics_figures/eGene_proportions.tiff", plot = p12, device = "tiff")

eQTL_eGene <- plot_grid(p11, p12, label_size = 1)

eQTL_eGene2 <- plot_grid(
  ggdraw() + draw_label("Proportion of eQTLs and eGenes in selection outliers", fontface = 'bold', size = 14, hjust = 0.5),
  eQTL_eGene,
  ncol = 1,
  rel_heights = c(0.1, 1)  # Adjust the height ratio as needed
)

ggsave("PLOS_genetics_figures/eQTL_and_eGenes_proportions.tiff", plot = eQTL_eGene2, device = "tiff", width = 10, height = 10)





# Change directory
setwd("../network_analysis/")

# plot connectivity
lfmm_connectivity <- read.table("lfmm_connectivity.txt")
names(lfmm_connectivity) <- c("gene", "connectivity")

pcadapt_connectivity <- read.table("pcadapt_connectivity.txt")
names(pcadapt_connectivity) <- c("gene", "connectivity")

eQTL_connectivity <- read.table("eQTL_connectivity.txt")
names(eQTL_connectivity) <- c("gene", "connectivity")

eGene_connectivity <- read.table("eGene_connectivity.txt")
names(eGene_connectivity) <- c("gene", "connectivity")

control_connectivity <- read.table("control_connectivity.txt")
names(control_connectivity) <- c("gene", "connectivity")

# Plot connectivity values for control and selection  outliers
eGene <- eGene_connectivity$connectivity
eQTL <- eQTL_connectivity$connectivity
LFMM <- lfmm_connectivity$connectivity
PCAdapt <- pcadapt_connectivity$connectivity
Control <- control_connectivity$connectivity
data <- cbind(eGene, eQTL, LFMM, PCAdapt, Control)

data <- data.frame(data) %>%
  mutate_all(as.numeric) %>%
  na.omit()


data_long <- data %>%
  gather(Method, connectivity) # %>%
#mutate(log_connectivity = log10(connectivity))

data_long <- data_long[is.finite(data_long$connectivity), ]

con_lfmm <- filter(data_long, Method == "LFMM")$connectivity
con_pcad <- filter(data_long, Method == "PCAdapt")$connectivity
con_control <- filter(data_long, Method == "Control")$connectivity
wilcox.test(con_control, con_lfmm) # (p-value < 2.2e-16, Wilcoxon rank sum test) 
wilcox.test(con_control, con_pcad) # ( p-value = 7.295e-05, Wilcoxon rank sum test) 
wilcox.test(con_pcad, con_lfmm) # p-value < 2.2e-16

# data3_long <- data3_long[is.finite(data3_long$FayH), ]

m <- median(filter(data_long, Method == "Control")$connectivity)

p13 <- data_long %>%
  mutate(connectivity = log(connectivity + 1))  %>% # Adding 1 to avoid log(0) issues
  ggplot(aes(x = Method, y = connectivity, fill = Method)) +
  geom_violin(adjust = 10, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     ref.group = "Control", label.y = 9) + # Specify label position
  #scale_y_continuous(trans = 'log10') +
  geom_hline(yintercept = 6.86, color = "black", alpha = 0.8) + 
  # scale_y_continuous(trans = 'log10') +# Position the horizontal line at log10(10000)
  xlab(" ") +
  ylab("In(connectivity)") +
  ggtitle(" ") +
  scale_fill_manual(values = my_colors) +
  guides(fill = "none") +
  theme_bw()


# Save the ggplot as a .tiff file
ggsave("../argo_plots/PLOS_genetics_figures/connectivity_plot.tiff", plot = p13, device = "tiff", width = 7, height = 4)


# Estimating population diversity stats
all_pop_theta2 <- all_pop_theta %>%
  dplyr::select(., tW, numSites) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  na.omit()
mean(all_pop_theta2$Watterson_theta)

all_pop_theta3 <- all_pop_theta %>%
  dplyr::select(., tP, numSites) %>%
  mutate(theta_pie = tP/numSites) %>%
  na.omit()
mean(all_pop_theta3$theta_pie)

mean(all_pop_theta$tajD)

all_pop_FST <- read.table(file="../../angsd/FST/fst_stat",sep="\t",header=F)
colnames(all_pop_FST) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")
mean(all_pop_FST$weighted_FST)

