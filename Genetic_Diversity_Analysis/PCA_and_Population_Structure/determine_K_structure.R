library(tidyverse)

# Function to load and format qopt files
read_qopt <- function(k) {
  q <- read.table(paste0("~/admixture/str_k_", k, ".qopt"))
  pop <- read.table("~/admixture/samples.txt")
  q <- data.frame(q, individuals = pop$V1)
  
  # Filter out bad samples
  exclude <- c("ARG1820.Aligned.sortedByCoord.out.bam", "ARG1805.Aligned.sortedByCoord.out.bam",
               "ARG1834.Aligned.sortedByCoord.out.bam", "Ames695.Aligned.sortedByCoord.out.bam",
               "Ames449.Aligned.sortedByCoord.out.bam", "btm7B-14.Aligned.sortedByCoord.out.bam",
               "btm5-1.Aligned.sortedByCoord.out.bam", "arg11B-11.Aligned.sortedByCoord.out.bam")
  q <- q %>% filter(!individuals %in% exclude)
  
  # Manually assign population label
  pop_labels <- c("Coast", "North", "North", "Coast", "North", "Coast", "Coast", "North", "Coast", "Coast",
                  "North", "North", "Coast", "Coast", "Coast", "North", "North", "North", "North")
  q$pop <- pop_labels
  
  # Reshape to long format for ggplot
  q_long <- q %>%
    arrange(pop) %>%
    mutate(ind = factor(individuals, levels = individuals)) %>%
    pivot_longer(cols = starts_with("V") | starts_with("K"),
                 names_to = "Cluster", values_to = "Proportion") %>%
    mutate(K = paste0("K = ", k))
  
  return(q_long)
}

# Read and combine data for K = 2, 3, 4
q2 <- read_qopt(2)
q3 <- read_qopt(3)
q4 <- read_qopt(4)
all_q <- bind_rows(q2, q3, q4)

# Custom color palettes
custom_colors <- list(
  "K = 2" = c("darkgrey", "brown"),
  "K = 3" = c("darkgrey", "brown", "aquamarine4"),
  "K = 4" = c("darkgrey", "brown", "aquamarine4", "gold")
)

plot_k <- function(q_data, k, colors, show_x_axis = FALSE) {
  q_k <- q_data %>% filter(K == paste0("K = ", k))
  
  ind_to_pop <- q_k %>%
    distinct(ind, pop) %>%
    arrange(ind) %>%
    deframe()
  
  ggplot(q_k, aes(x = ind, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black", size = 0.2) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = ind_to_pop) +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x = if (show_x_axis) element_text(angle = 90, vjust = 0.5, hjust = 1, size = 22) else element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_text(size = 22),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = if (show_x_axis) margin(5, 5, 35, 55) else margin(5, 5, 5, 55)  # ← Left = 60, Bottom = 40 only for K=4
    ) +
    ylab(paste0("K = ", k))
}

p2 <- plot_k(all_q, 2, custom_colors[["K = 2"]])
p3 <- plot_k(all_q, 3, custom_colors[["K = 3"]])
p4 <- plot_k(all_q, 4, custom_colors[["K = 4"]], show_x_axis = TRUE)

(p2 / p3 / p4) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(0, 0, 0, 0), plot.spacing = unit(0, "pt"))




#-------------------------------------------------------------------------------
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

write.table(covariates, file = "eQTL_analysis/covariates.txt", sep = "\t", quote = F)
