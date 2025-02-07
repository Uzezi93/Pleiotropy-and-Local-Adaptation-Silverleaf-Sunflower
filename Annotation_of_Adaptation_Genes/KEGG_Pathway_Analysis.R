setwd("~/Documents/sunflower/")
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
# library(AnnotationHub) 
# library(GenomicFeatures)
library(rtracklayer)
library(pRoloc)
library(pathview)

#KEGG overreprentation of gene modules with more than expected lfmm and pcadapt genes
search_kegg_organism('han', by='kegg_code')

han <- search_kegg_organism('Helianthus annuus', by='scientific_name')
dim(han)

#------------------------------------------------------------
# Get ncbi gene IDs for LFMM genes
lfmm_gff <- "lfmm_exons.gff"

# Read the GFF file
lfmm_gff <- readGFF(lfmm_gff)

# Extract gene IDs
lfmm_gff_genes <- lfmm_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = lfmm_gff_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q1 <- barplot(kk) +
  ggtitle("KEGG Pathway Over Representation Analysis for LFMM Genes")

ggsave("LFMM_Pathway_KEGG_pathway_ORA_barplot.tiff", plot = q1, device = "tiff")

q2 <- dotplot(kk) +
  ggtitle("LFMM Module KEGG Pathway ORA")

ggsave("LFMM_Pathway_KEGG_pathway_OVA_dotplot.tiff", plot = q2, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_LFMM_attribute.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# Get ncbi gene IDs PCAdapt outliers
pcadapt_gff <- "pcadapt_exons.gff"

# Read the GFF file
pcadapt_gff <- readGFF(pcadapt_gff)

# Extract gene IDs
pcadapt_gff_genes <- pcadapt_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = pcadapt_gff_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q3 <- barplot(kk) +
  ggtitle("KEGG Pathway Over Representation Analysis for PCAdapt Genes")

ggsave("PCAdapt_Pathway_KEGG_pathway_ORA_barplot.tiff", plot = q3, device = "tiff")

q4 <- dotplot(kk) +
  ggtitle("PCAdapt Module KEGG Pathway ORA")

ggsave("PCAdapt_Pathway_KEGG_pathway_OVA_dotplot.tiff", plot = q4, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_PCAdapt_attribute.txt", sep = "\t", quote = F)

#-----------------------------------------------------------
# Get ncbi gene IDs for genes identified by both PCAdapt and LFMM methods
lfmm_pcadapt_gff <- "lfmm_pcadapt_exons.gff"

# Read the GFF file
lfmm_pcadapt_gff <- readGFF(lfmm_pcadapt_gff)

# Extract gene IDs
lfmm_pcadapt_gff_genes <- lfmm_pcadapt_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = lfmm_pcadapt_gff_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q5 <- barplot(kk) +
  ggtitle("KEGG Pathway ORA For Genes Identified by LFMM and PCAdapt")

ggsave("LFMM_PCAdapt_Pathway_KEGG_pathway_OVA.tiff", plot = q5, device = "tiff")

q6 <- dotplot(kk) +
  ggtitle("KEGG Pathway ORA For Genes Identified by LFMM and PCAdapt")

ggsave("LFMM_PCAdapt_Pathway_KEGG_pathway_OVA_dotplot.tiff", plot = q6, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_LFMM_PCAdapt_attribute.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# Get ncbi gene IDs for LFMM Darkgrey Module
darkgrey_gff <- "LFMM_darkgrey_exons.gff"

# Read the GFF file
darkgrey_gff <- readGFF(darkgrey_gff)

# Extract gene IDs
darkgrey_genes <- darkgrey_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = darkgrey_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q7 <- barplot(kk) +
  ggtitle("Darkgrey Module KEGG Pathway ORA")

ggsave("Darkgrey_Pathway_KEGG_pathway_OVA.tiff", plot = q7, device = "tiff")

q8 <- dotplot(kk) +
  ggtitle("Darkgrey Module KEGG Pathway ORA")

ggsave("Darkgrey_KEGG_pathway_OVA_dotplot.tiff", plot = q8, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_Darkgrey_attribute.txt", sep = "\t", quote = F)

#-----------------------------------------------------------
# KEGG Module over representation analysis darkgrey module
kk <- enrichMKEGG(gene         = darkgrey_genes,
                  organism     = 'han',
                  pvalueCutoff = 0.05)
head(kk) 

y1 <- barplot(kk) +
  ggtitle("Darkgrey KEGG Module ORA")

module <- data.frame(kk)

write.table(module, file = "KEGG_module_ORA_Darkgrey.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# Get ncbi gene IDs for LFMM Darkturqoiuse Module
darkturqouise_gff <- "LFMM_darkturquoise_exons.gff"

# Read the GFF file
darkturqouise_gff <- readGFF(darkturqouise_gff)

# Extract gene IDs
darkturqouise_genes <- darkturqouise_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = darkturqouise_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q9 <- barplot(kk) +
  ggtitle("Darkturqouise Module KEGG Pathway ORA")

ggsave("Darkturqouise_Pathway_KEGG_pathway_OVA.tiff", plot = q9, device = "tiff")

q10 <- dotplot(kk) +
  ggtitle("Darkturqouise Module KEGG Pathway ORA")

ggsave("Darkturqouise_KEGG_pathway_OVA_dotplot.tiff", plot = q10, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_Darkturqouise_attribute.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# KEGG Module over representation analysis darkturqouise module
kk <- enrichMKEGG(gene         = darkturqouise_genes,
                  organism     = 'han',
                  pvalueCutoff = 0.05)
head(kk) 

y2 <- barplot(kk) +
  ggtitle("Darkturqouise KEGG Module ORA")

module <- data.frame(kk)

write.table(module, file = "KEGG_module_ORA_Darkturqouise.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# Get ncbi gene IDs for LFMM Green Module
green_gff <- "LFMM_green_exons.gff"

# Read the GFF file
green_gff <- readGFF(green_gff)

# Extract gene IDs
green_genes <- green_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = green_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

# Repair the names of pathway description
kk@result$Description <- gsub(" - Helianthus annuus \\(common sunflower\\)", "", kk@result$Description)

head(kk@result$Description)

q11 <- barplot(kk) +
  ggtitle("Green Module KEGG Pathway ORA")

ggsave("Green_Pathway_KEGG_pathway_OVA.tiff", plot = q11, device = "tiff")

q12 <- dotplot(kk) +
  ggtitle("Green Module KEGG Pathway ORA")

ggsave("Green_KEGG_pathway_OVA_dotplot.tiff", plot = q12, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "KEGG_pathway_ORA_Green_attribute.txt", sep = "\t", quote = F)

#------------------------------------------------------------
# KEGG Module over representation analysis darkturqouise module
kk <- enrichMKEGG(gene         = green_genes,
                  organism     = 'han',
                  pvalueCutoff = 0.05)
head(kk) 

y3 <- barplot(kk) +
  ggtitle("Green KEGG Module ORA")

module <- data.frame(kk)

write.table(module, file = "KEGG_module_ORA_Green.txt", sep = "\t", quote = F)




