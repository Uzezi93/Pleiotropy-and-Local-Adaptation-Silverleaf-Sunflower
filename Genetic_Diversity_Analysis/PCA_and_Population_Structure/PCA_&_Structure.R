setwd("Argo_redo/")

library(tidyverse)
library(caret)
library(clusterSim)
library(FactoMineR)
library(factoextra)
library(pls)

# read in PCA data
pca <- read_table("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/argo.eigenvec", col_names = FALSE)
eigenval <- scan("/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/argo.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

coords <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/arg_native_poplatlong.csv")

# Sample coordinates (latitude)
coord2 <- coords %>%
  filter(pop %in% pop[1:22]) %>%
  filter(!pop %in% c("arg11B", "btm5", "btm7B"))


# Step 1: Fix column types in pca (they're character now)
pca <- pca %>%
  filter(ind != "IID") %>%            # remove header row
  mutate(across(starts_with("PC"), as.numeric))

# Step 2: Harmonize IDs between datasets
# coords$pop has IDs like "btm9", but pca2$ind has "btm9-4".
# Let's strip replicate suffixes (like -4) to match coords.
pca <- pca %>%
  mutate(base_id = sub("-[0-9]+$", "", ind))

# Step 3: Join with coords
pca <- pca %>%
  inner_join(coords, by = c("base_id" = "pop"))

pca <- pca %>%
  rename(population = X)

# remake data.frame
row.names(pca) <- pca$ind
pca <- pca[order(row.names(pca)),]
n <- row.names(pca)
pca$ind <- NULL

write.table(pca, file = "~/genetic_pca.txt", sep = "\t", quote = F, row.names = T)

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
    axis.title.x = element_text(size = 20),   # X-axis label
    axis.title.y = element_text(size = 20),   # Y-axis label
    axis.text.x  = element_text(size = 20),   # X-axis ticks
    axis.text.y  = element_text(size = 20),   # Y-axis ticks
    legend.title = element_text(size = 20),   # Legend title
    legend.text  = element_text(size = 20),   # Legend labels
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5)  # Optional plot title
  )


ggsave(
  plot = c,
  filename = "/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision/PCA.tiff",
  width = 7,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"  # lossless compression required by many journals
)

# Top two samples from PC2
pca2 <- pca %>% 
  rownames_to_column(var = "ind")   # turn rownames into a real column

pca2 %>%
  dplyr::select(ind, PC1, PC2, population, base_id) %>%
  arrange(desc(PC2)) %>%
  slice(1:2)

library(tidyverse)
library(stringr)

admix_dir <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/NGSADMIX"
fam_path  <- "/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/argo_snps.fam"

# ------------ load/prepare metadata ------------
coord <- read.csv("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/arg_native_poplatlong.csv",
                  header = TRUE) %>%
  filter(X %in% c("Coast","North")) %>%
  transmute(sample = pop, group = X)                      # sample = base id; group = Coast/North

fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE) |>
  setNames(c("FID","IID","PID","MID","SEX","PHENO")) |>
  mutate(sample = sub("-.*", "", IID))                    # strip suffix to base id

# join to attach group to each IID
meta <- fam %>%
  left_join(coord, by = "sample")

# check for missing group assignments
if (any(is.na(meta$group))) {
  warning("Some fam samples missing in coord: ",
          paste(unique(meta$sample[is.na(meta$group)]), collapse = ", "))
}

# build the desired order:
#   1) Coast before North
#   2) within each, follow the order in coord$sample (left-to-right in your coord file)
coord_order <- coord$sample
meta <- meta %>%
  mutate(group = factor(group, levels = c("Coast","North")),
         sample_ord = factor(sample, levels = coord_order))

ord <- order(meta$group, meta$sample_ord, na.last = TRUE)

# vectors in final plotting order
ids   <- meta$IID[ord]
pops  <- meta$sample[ord]   # base sample names used for boundaries/labels
group <- as.character(meta$group[ord])

# ------------ read all Q files & reshape ------------
qfiles <- list.files(admix_dir, pattern = "\\.admix\\.[0-9]+\\.Q$", full.names = TRUE)
stopifnot(length(qfiles) > 0)
getK  <- function(f) as.integer(stringr::str_match(basename(f), "\\.admix\\.(\\d+)\\.Q$")[,2])

df <- purrr::map_dfr(qfiles, function(f){
  K <- getK(f)
  Q <- as.matrix(read.table(f, header = FALSE, stringsAsFactors = FALSE))
  Q <- Q[ord, , drop = FALSE]                              # reorder rows to our meta order
  colnames(Q) <- paste0("C", seq_len(ncol(Q)))
  as_tibble(Q) %>%
    mutate(IID = ids, pop = pops, group = group, K = K) %>%
    pivot_longer(starts_with("C"), names_to = "cluster", values_to = "q")
}) %>%
  mutate(
    IID = factor(IID, levels = ids),                       # lock individual order
    K   = factor(K, levels = sort(unique(K)), labels = paste0("K=", sort(unique(K))))
  )

df_sub <- df %>% filter(K %in% c("K=2","K=3","K=4"))

# ------------ population boundaries & headers ------------
# First index of each pop (base id) along the x-axis of IIDs
pop_boundaries <- tibble(IID = ids, pop = pops) %>%
  group_by(pop) %>% summarise(xpos = which(ids == first(IID))[1], .groups = "drop")

# where to draw Coast/North bars above
n_coast <- sum(group == "Coast")
n_total <- length(ids)
pop_groups <- tibble(
  label  = c("Coast","North"),
  xstart = c(1, n_coast + 1),
  xend   = c(n_coast, n_total),
  y      = 1.05
)

# === build label positions (add clean label without "-suffix") ===
iid_pos <- tibble(
  IID   = ids,
  xpos  = seq_along(ids),
  label = sub("-.*", "", ids)   # <- remove "-" and everything after
)

# === plotting ===
# make the custom labels bigger and slightly lower
label_size <- 4.8

p <- ggplot(df_sub, aes(x = IID, y = q, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., scales = "free_x", space = "free_x") +
  geom_segment(data = pop_boundaries,
               aes(x = xpos - 0.5, xend = xpos - 0.5, y = 0, yend = 1),
               inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  # x-axis labels only for K=4 (bigger)
  geom_text(
    data = df_sub %>%
      filter(K == "K=4") %>%
      distinct(IID, K) %>%
      left_join(iid_pos, by = "IID"),
    aes(x = xpos, y = -0.09, label = label),   # lower a bit
    inherit.aes = FALSE, angle = 90, vjust = 0.5, size = label_size
  ) +
  geom_segment(data = pop_groups,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = pop_groups,
            aes(x = (xstart + xend)/2, y = y + 0.1, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 4) +
  labs(x = "Individuals (ordered by population)", y = "Ancestry proportion") +
  scale_fill_manual(values = c("brown","aquamarine4","gold","grey")) +
  coord_cartesian(ylim = c(-0.12, 1.2), clip = "off") +   # more room for labels
  theme_minimal(base_size = 16) +
  theme(
    legend.position   = "none",
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.title.x      = element_text(size = 16, face = "bold", margin = margin(t = 12)),
    axis.title.y      = element_text(size = 16, face = "bold"),
    axis.text.y       = element_text(size = 14),
    panel.grid.major.x= element_blank(),
    panel.grid.minor  = element_blank(),
    strip.text.y      = element_text(face = "bold", size = 12),
    plot.margin       = margin(5.5, 5.5, 65, 5.5)  # extra bottom space
  )

print(p)


grid_fig <- cowplot::plot_grid(
  c, p,
  labels = c("A", "B"),
  label_size = 20,   # increase label size (default is 14)
  ncol = 2,
  align = "hv"
)


ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/argo_revision/pca_structure.tiff", grid_fig, width = 20, height = 10, dpi = 300)
