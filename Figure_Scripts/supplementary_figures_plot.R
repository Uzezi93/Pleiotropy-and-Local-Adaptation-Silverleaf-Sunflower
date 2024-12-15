setwd("../argo_plots/")

# Plot supplementary figures
# 1. Relationship between outlier pvalues and connectivity
# write.table(argo.pv, file = "./argo.pv.txt", sep = "\t", quote = FALSE)
# write.table(argo.pv2, file = "argo.pv2.txt", sep = "\t", quote = FALSE)
# write.table(argo.pv3, file = "argo.pv3.txt", sep = "\t", quote = FALSE)
# write.table(argo.pv4, file = "argo.pv4.txt", sep = "\t", quote = FALSE)
# write.table(argo.pv5, file = "argo.pv5.txt", sep = "\t", quote = FALSE)
# write.table(res, file = "pcadapt_res.txt", sep = "\t", quote = FALSE)


#LFMM_pvalue <- data.frame(cbind(argo.pv$pvalue, argo.pv2$pvalue, argo.pv3$pvalue, argo.pv4$pvalue, argo.pv5$pvalue)) %>%
 # na.omit()
#LFMM_pvalue <- gather(LFMM_pvalue)
#LFMM_pvalue <- LFMM_pvalue$value
#PCAdapt_pvalue <- res$pvalues

# Save to pvalues
#cat(PCAdapt_pvalue, file = "PCAdapt_pvalues.txt")
#cat(LFMM_pvalue, file = "LFMM_pvalues.txt")

# Get connectivity values
# Change directory
LFMM_pvalue <- read.table("../network_analysis/lfmm_pvalues.txt", header = F)
names(LFMM_pvalue) <- c("pvalue", "gene")
PCAdapt_pvalue <- read.table("../network_analysis/pcadapt_pvalues.txt", header = F)
names(PCAdapt_pvalue) <- c("pvalue", "gene")

setwd("../network_analysis")

# plot connectivity
lfmm_connectivity <- read.table("connectivity.txt", header = T)
names(lfmm_connectivity) <- c("focal_gene", "connectivity")

lfmm_connectivity <- lfmm_connectivity %>%
  separate(focal_gene, into = c("name", "gene")) %>%
  select(gene, connectivity) %>%
  distinct()

# Merge by "gene"
lfmm <- merge(LFMM_pvalue, lfmm_connectivity, by = c("gene")) 
lfmm$connectivity <- as.numeric(lfmm$connectivity)

lfmm <- lfmm %>%
  mutate(log_connectivity = log(connectivity))

pcadapt_connectivity <- read.table("connectivity.txt", header = T)
names(pcadapt_connectivity) <- c("focal_gene", "connectivity")

pcadapt_connectivity <- pcadapt_connectivity %>%
  separate(focal_gene, into = c("name", "gene")) %>%
  select(gene, connectivity) %>%
  distinct()

# Merge by "gene"
pcadapt <- merge(PCAdapt_pvalue, pcadapt_connectivity, by = c("gene")) 
pcadapt$connectivity <- as.numeric(pcadapt$connectivity)

pcadapt <- pcadapt %>%
  mutate(log_connectivity = log(connectivity))

# Perform linear regression
model <- lm(pvalue ~ log_connectivity, data = lfmm)

# Display summary of the model
summary(model)

# Create scatter plot with regression line
lfmm_lm_connectivity <- ggplot(lfmm, aes(x = log_connectivity, y = pvalue)) +
  #geom_point() +
  geom_smooth(method = "gam", col = "red") +
  stat_cor(method = "pearson", label.y = 0.53) +
  labs(x = "Connectivity", y = "LFMM p-values") +
  theme_bw()


# Perform linear regression
model <- lm(pvalue ~ log_connectivity, data = pcadapt)

# Display summary of the model
summary(model)

# Create scatter plot with regression line
pcadapt_lm_connectivity <- ggplot(pcadapt, aes(x = log_connectivity, y = pvalue)) +
  #geom_point() +
  geom_smooth(method = "gam", col = "red") +
  stat_cor(method = "pearson", label.y = 0.55) +
  labs(x = "Connectivity", y = "PCAdapt p-values") +
  theme_bw()

lm_connectivity <- plot_grid(lfmm_lm_connectivity, pcadapt_lm_connectivity)
ggsave("PLOS_genetics_figures/pval_lm_connectivity.tiff", plot = lm_connectivity, device = "tiff", width = 10, height = 5)


# density plot of connectivity values
lfmm_connectivity_density <- ggplot(data, aes(x = lfmm_connectivity)) +
  geom_density() +
  labs(x = "Connectivity", y = "LFMM p-values")

pcadapt_connectivity_density <- ggplot(data2, aes(x = pcadapt_connectivity)) +
  geom_density() +
  labs(x = "Connectivity", y = "PCAdapt p-values")

density_connectivity <- plot_grid(lfmm_connectivity_density , pcadapt_connectivity_density)
ggsave("PLOS_genetics_figures/connectivity_density.tiff", plot = density_connectivity, device = "tiff", width = 10, height = 5)


# Plot pairwise fst
fst <- read.table(file="../../angsd/FST/fst_stat",sep="\t",header=F);
colnames(fst) <- c("chromosomes", "positions", "unweighted_FST", "weighted_FST")

# Step 1: Load the data
filename <- "../../angsd/FST/coastal.north_inland.ml"

##run in R
cn <-scan("../../angsd/FST/coastal.north_inland.ml")
source("plot2dSFS.R")
plot2<-function(s,...){
  dim(s)<-c(21,19)
  s[1]<-NA
  s[21,19]<-NA
  s<-s/sum(s,na.rm=T)
  
  pal <- color.palette(c("aquamarine4","brown"), space="rgb")
  pplot(s/sum(s,na.rm=T),pal=pal,...)
}

plot2(cn,ylab="coast",xlab="north")

fst <- read.table("../../angsd/FST/slidingwindow", skip=1, header=FALSE)
names(fst)[c(3,5)] <- c("midp", "fst")
plot(fst$midp, fst$fst, xlab="Physical position", ylab="Fst", col="#5f9ea0", pch=16)







# Do angsd for north and coast populations to get pestg files
north_pop_theta <- read.table(file = "../../angsd/SFS2/north_pop.thetas.idx.pestPG", sep = "\t", header = F)

coast_pop_theta <- read.table(file = "../../angsd/SFS2/coast_pop.thetas.idx.pestPG", sep = "\t", header = F)

# Rename columns
colnames(north_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")
colnames(coast_pop_theta) <- c("(indexStart,indexStop)", "chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

north_pop <- north_pop_theta %>%
  dplyr::select(., tW, numSites, tP, tajD) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  mutate(Theta_pie = tP/numSites) %>%
  na.omit()
tajD <- mean(north_pop$tajD) # -0.1745136
Watterson_theta <- mean(north_pop$Watterson_theta) # 0.2076168
Theta_pie <- mean(north_pop$Theta_pie) # 0.1937444

coast_pop <- coast_pop_theta %>%
  dplyr::select(., tW, numSites, tP, tajD) %>%
  mutate(Watterson_theta = tW/numSites) %>%
  mutate(Theta_pie = tP/numSites) %>%
  na.omit()
tajD <- mean(coast_pop$tajD) # -0.2985429
Watterson_theta <- mean(coast_pop$Watterson_theta) # 0.2353111
Theta_pie <- mean(coast_pop$Theta_pie) # 0.2060006




