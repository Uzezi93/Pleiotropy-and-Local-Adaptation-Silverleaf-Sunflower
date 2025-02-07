library(BiocManager)
library(reshape)
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(GenomicFeatures)
library(vsn)
library(EDASeq)
library(affy)
library(gplots)
library(ggpmisc)


# Load data
counts <- read.table("counts.txt",sep="\t",header=TRUE,stringsAsFactors=F)

# Examine dataset
str(counts)

# Convert count data to integers.
counts[, -c(1:6)] <- lapply(counts[, -c(1:6)], as.integer)

# Remove duplicate rows
counts <- counts %>% dplyr::distinct()

# backup original names
orig_names <- names(counts) # keep a back-up copy of the original names

# Change names
names(counts)[7:33] <- c("Ames449", "Ames695", "arg11B-11", "arg14B-7", "ARG1805", "ARG1820", "ARG1834", "arg2B-4", "arg4B-8", "arg6B-1", "btm10-5", "btm13-4", "btm17-4", "btm19-1", "btm20-8", "btm21-4", "btm22-8", "btm25-2", "btm26-4", "btm27-3", "btm30-6", "btm31-6", "btm32-3", "btm34-6", "btm5-1", "btm7B-14", "btm9-4")

# Examine new manipulated data
str(counts)

# Remove genebank accessions
counts$Ames449 <- NULL
counts$Ames695 <- NULL
counts$ARG1805 <- NULL
counts$ARG1820 <- NULL
counts$ARG1834 <- NULL

#DeSeq2 setup
## gene IDs should be stored as row.names
row.names(counts) <- make.names(counts$Geneid, unique = TRUE)

# Remove genes without counts and save loci for downstream analysis.
gene_loc <- counts %>%
  dplyr::select(c(1:4))

# write.table(gene_loc, quote = FALSE, file = "Helianthus_argophyllus project/gene_loc.txt")

## exclude the columns without read counts (columns 1 to 6 contain additional
## info such as genomic coordinates)
counts <- counts[ , -c(1:6)]

head(counts)

str(counts)

# let's use the info from our readcounts object

# Populations for score plot
# 1.	Coastal (samples arg6B-1, btm9-4, btm10-5, btm13-4, btm17-4, btm19-1, btm 20-8, btm21-4, btm22-8)

# 2. North Inland (samples arg14B-7, arg2B-4, arg4B-8, btm25-2, btm26-4, btm27-3, btm30-6, btm31-6, btm32-3, and btm34-6)

# 3. South Inland (samples arg11B-11, btm5-1, btm7B-14
a <- names(counts)

condition <- c("South Inland", "North Inland", "North Inland", "North Inland", "Coastal", "Coastal", "Coastal", "Coastal", "Coastal", "Coastal", "Coastal", "Coastal", "North Inland", "North Inland", "North Inland", "North Inland", "North Inland", "North Inland", "North Inland", "South Inland", "South Inland", "Coastal")

coldata <- data.frame(cbind(a,condition))
coldata$condition <- factor(coldata$condition)
row.names(coldata) <- coldata$a
coldata$a <- NULL
coldata$condition <- gsub(" ", "_",coldata$condition)


# sample_info <- DataFrame(condition = gsub("_[0-9]+", "", names(counts)),
#                         row.names = names(counts))
# sample_info

# str(sample_info)


# Create a coldata frame: its rows correspond to columns of data (i.e., matrix representing the countData)
#coldata <- data.frame(row.names=colnames(counts), sample)

# Generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition)
DESeq.ds

head(counts(DESeq.ds))

#------------------------------------------------------------------------
x <- counts(DESeq.ds)
x <- as.data.frame(x)
str(x)
head(x)
x <- melt(x)
my.formula <- y ~ x
p <- ggplot(data = x, aes(x = variable, y = value)) +
geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point()

p

#How many reads were sequenced for each sample ( = library sizes)?
colSums(counts(DESeq.ds))

colSums(counts(DESeq.ds)) %>% barplot

#------------------------------------------------------------------------

# Remove genes with no reads.
dim(DESeq.ds)

keep_genes <- rowSums(counts(DESeq.ds)) > 0
DESeq.ds <- DESeq.ds[ keep_genes, ]
dim(DESeq.ds)

# As you can see, there are now fewer features stored in the DESeq.ds (first entry of the dim() result). 

# The filtering was also translated to the count matrix that we store in that object (and all other matrices stored in the assay slot).

counts(DESeq.ds) %>% str

assay(DESeq.ds) %>% str

#------------------------------------------------------------------------

# Normalizing for sequencing depth and RNA composition difference

## define a function to calculate the geometric mean
gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

## calculate the geometric mean for each gene using that function
## note the use of apply(), which we instruct to apply the gm_mean()
## function per row (this is what the second parameter, 1, indicates)
pseudo_refs <- counts(DESeq.ds) %>% apply(., 1, gm_mean)

## divide each value by its corresponding pseudo-reference value
pseudo_ref_ratios <- counts(DESeq.ds) %>% apply(., 2, function(cts){ cts/pseudo_refs})

## if you want to see what that means at the single-gene level,
## compare the result of this:
counts(DESeq.ds)[1,]/pseudo_refs[1]

## with
pseudo_ref_ratios[1,]

## determine the median value per sample to get the size factor
apply(pseudo_ref_ratios , 2, median)


# Calculating and applying the size factor
DESeq.ds <- estimateSizeFactors(DESeq.ds) # calculate SFs, add them to object
plot( sizeFactors(DESeq.ds), colSums(counts(DESeq.ds)), # assess them
      ylab = "library sizes", xlab = "size factors", cex = .6 )

## setting up the plotting layout
par(mfrow=c(1,2))
## extracting normalized counts
counts.sf_normalized <- counts(DESeq.ds, normalized=TRUE)
## adding the boxplots
boxplot(counts.sf_normalized, main = "SF normalized", cex = .6)
boxplot(counts(DESeq.ds), main = "read counts only", cex = .6)

# To see the influence of the sequencing depth normalization, make two box plots of log2(read counts):
#- one for non-normalized counts - the other one for normalized counts

par(mfrow=c(1,2)) # to plot the two box plots next to each other
## bp of non-normalized
boxplot(log2(counts(DESeq.ds)+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts)", cex = .6)
## bp of size-factor normalized values
boxplot(log2(counts(DESeq.ds, normalize= TRUE) +1), notch=TRUE,
        main = "Size-factor-normalized read counts",
        ylab="log2(read counts)", cex = .6)

#------------------------------------------------------------------------

library(reshape2)
library(ggpmisc)
a <- log2(counts(DESeq.ds)+1)
a <- data.frame(a)
a = melt(a)
head(a)

dev.off()

my.formula <- y ~ x
p <- ggplot(data = a, aes(x = variable, y = value)) +
  geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point()

print(p)

## generate the base meanSdPlot using sequencing depth normalized log2(read counts)
log.norm.counts <- log2(counts(DESeq.ds, normalized=TRUE) + 1)
## set up ploting frames
par(mfrow=c(1,1))
## generate the plot
msd_plot <- vsn::meanSdPlot(log.norm.counts,
                            ranks=FALSE, # show the data on the original scale
                            plot = FALSE)
## since vsn::meanSdPlot generates a ggplot2 object, this can be
## manipulated in the usual ways
msd_plot$gg +
  ggtitle("Sequencing depth normalized log2(read counts)") +
  ylab("standard deviation")

#------------------------------------------------------------------------
# Reducing the dependence of the variance on the mean
## this actually generates a different type of object!
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
## set blind = FALSE if the conditions are expected to introduce
## strong differences in a large proportion of the genes

# Let's visually check the results of the rlog transformation:
par(mfrow=c(1,2))
plot(log.norm.counts[,1:2], cex=.1,
     main = "size factor and log2-transformed")

## the rlog-transformed counts are stored in the accessor "assay"
plot(assay(DESeq.rlog)[,1],
     assay(DESeq.rlog)[,2],
     cex=.1, main = "rlog transformed",
     xlab = colnames(assay(DESeq.rlog[,1])),
     ylab = colnames(assay(DESeq.rlog[,2])) )

rlog.norm.counts <- assay(DESeq.rlog)
rlog.norm.counts <- data.frame(rlog.norm.counts)
write.table(rlog.norm.counts, quote = F, file="Helianthus_argophyllus project/normalized_readcounts.txt")

non_norm <- assay(DESeq.ds)
# What does the mean-sd-plot show?
## rlog-transformed read counts
msd_plot <- vsn::meanSdPlot(Srlog.norm.counts, ranks=FALSE, plot = FALSE)
msd_plot$gg + ggtitle("rlog transformation") + coord_cartesian(ylim = c(0,3))

##-----------------------------------------------------------------------

# Regularized log transformation for clustering/heatmaps, etc
head(assay(DESeq.rlog))

# Compare normalized and non-normalized data
par(mfrow=c(1,2))
hist(assay(DESeq.rlog))
hist(assay(DESeq.ds))

par(mfrow=c(1,1))

# Rename samples according to their populations
# Populations for score plot
# 1.	Coastal (samples arg6B-1, btm9-4, btm10-5, btm13-4, btm17-4, btm19-1, btm 20-8, btm21-4, btm22-8)

# 2. North Inland (samples arg14B-7, arg2B-4, arg4B-8, btm25-2, btm26-4, btm27-3, btm30-6, btm31-6, btm32-3, and btm34-6)

# 3. South Inland (samples arg11B-11, btm5-1, btm7B-14
# Principal Components Analysis for grouped samples

plotPCA(DESeq.rlog)

# Colors for plots below
## Ugly:
(mycols <- 1:length(unique(condition)))

## Use RColorBrewer, better
library(RColorBrewer)

mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(DESeq.rlog))))
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=T, trace="none",
          col=colorpanel(100, "black", "green", "red"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")

dev.off()

# Get differential expression results
# res <- results(DESeq(DESeq.ds))
# table(res$padj<0.05)

#------------------------------------------------------------------------















