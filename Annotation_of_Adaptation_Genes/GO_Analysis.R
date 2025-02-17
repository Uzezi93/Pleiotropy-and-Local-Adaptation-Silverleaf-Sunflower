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

#------------------------------------------------------------
# Get ncbi gene IDs for LFMM and PCAdapt outliers
lfmm_gff <- "lfmm_exons.gff"

# Read the GFF file
lfmm_gff  <- readGFF(lfmm_gff )

# Extract gene IDs
lfmm_genes <- lfmm_gff$gene 

#------------------------------------------------------------

pcadapt_gff <- "pcadapt_exons.gff"

# Read the GFF file
pcadapt_gff <- readGFF(pcadapt_gff)

# Extract gene IDs
pcadapt_genes <- pcadapt_gff$gene

#------------------------------------------------------------

# Get ncbi gene IDs for genes identified by LFMM and PCAdapt
lfmm_pcadapt_gff <- "lfmm_pcadapt_exons.gff"

# Read the GFF file
lfmm_pcadapt_gff  <- readGFF(lfmm_pcadapt_gff)

# Extract gene IDs
lfmm_pcadapt_genes <- lfmm_pcadapt_gff$gene 

#-----------------------------------------------------------
# Get ncbi gene IDs for genes in the LFMM darkgrey module
darkgrey_gff <- "LFMM_darkgrey_exons.gff"

# Read the GFF file
darkgrey_gff  <- readGFF(darkgrey_gff)

# Extract gene IDs
darkgrey_genes <- darkgrey_gff$gene 

#-----------------------------------------------------------
# Get ncbi gene IDs for genes in the LFMM darkturqoise module
darkturqouise_gff <- "LFMM_darkturquoise_exons.gff"

# Read the GFF file
darkturqouise_gff  <- readGFF(darkturqouise_gff)

# Extract gene IDs
darkturqouise_genes <- darkturqouise_gff$gene 

#------------------------------------------------------------
# Get ncbi gene IDs for genes in the LFMM green module
green_gff <- "LFMM_green_exons.gff"

# Read the GFF file
green_gff  <- readGFF(green_gff)

# Extract gene IDs
green_genes <- green_gff$gene 
#------------------------------------------------------------
# Get ncbi gene IDs for genes in the eGenes
eGene_gff <- "eGenes_exons.gff"

# Read the GFF file
eGene_gff  <- readGFF(eGene_gff)

# Extract gene IDs
eGene_genes <- eGene_gff$gene 
#------------------------------------------------------------

all_genes <- read.delim("~/Downloads/GCF_002127325.2_HanXRQr2.0-SUNRISE_gene_ontology.gaf", skip = 8, header = TRUE)

# filter all_genes for genes with non-zero read counts


# Get GO description for all genes
go_terms <- clusterProfiler::go2term(all_genes$GO_ID)

go_terms <- go_terms %>% 
  dplyr::rename(GO_ID = go_id)

all_genes <- left_join(all_genes, go_terms, by = "GO_ID")

term2gene=all_genes[, c("Term", "Symbol")] %>%
  na.omit

x = enricher(lfmm_genes, TERM2GENE=term2gene)
head(summary(x))

# Save gene over representation analysis table
write.table(x, file = "LFMM_Gene_over-representation.txt", sep = "\t", quote = F)

p1 <- barplot(x, main="LFMM Genes Over Representation Analysis") +
  ggtitle("Gene Over Representation Analysis for LFMM Genes")

# Save the ggplot as a .tiff file
ggsave("LFMM_gene_over_representation_barplot.tiff", plot = p1, device = "tiff", width = 7, height = 4)

p2 <- dotplot(x) +
  ggtitle("Gene Over Representation Analysis for LFMM Genes")
ggsave("LFMM_gene_over_representation_dotplot.tiff", plot = p2, device = "tiff", width = 7, height = 4)

#------------------------------------------------------------
# Enrichnment for PCAdapt Genes
y = enricher(pcadapt_genes, TERM2GENE=term2gene)
head(summary(y))

# Save gene over representation analysis table
write.table(y, file = "PCAdapt_Gene_over-representation.txt", sep = "\t", quote = F)

p3 <- barplot(y, main="PCAdapt Genes Over Representation Analysis") +
  ggtitle("Gene Over Representation Analysis for PCAdapt Genes")

# Save the ggplot as a .tiff file
ggsave("PCAdapt_gene_over_representation_barplot.tiff", plot = p3, device = "tiff", width = 7, height = 4)

p4 <- dotplot(y) +
  ggtitle("Gene Over Representation Analysis for PCAdapt Genes")
ggsave("PCAdapt_gene_over_representation_dotplot.tiff", plot = p4, device = "tiff", width = 7, height = 4)

#------------------------------------------------------------
# Enrichment analysis for genes identified by LFMM and PCAdapt
# Enrichnment for PCAdapt Genes
z = enricher(lfmm_pcadapt_genes, TERM2GENE=term2gene)
head(summary(z))

# Save gene over representation analysis table
write.table(z, file = "LFMM_PCAdapt_Gene_over-representation.txt", sep = "\t", quote = F)

p5 <- barplot(z, main="Over Representation Analysis for LFMM and PCAdapt genes") +
  ggtitle("Gene ORA for Genes Identied by LFMM and PCAdapt")

# Save the ggplot as a .tiff file
ggsave("LFMM_PCAdapt_gene_over_representation_barplot.tiff", plot = p5, device = "tiff", width = 7, height = 4)

p6 <- dotplot(z) +
  ggtitle("Gene Over Representation Analysis for Genes identied by LFMM and PCAdapt Methods")
ggsave("LFMM_PCAdapt_gene_over_representation_dotplot.tiff", plot = p6, device = "tiff", width = 7, height = 4)

#------------------------------------------------------------
# Enrichment analysis for genes LFMM darkgrey module
a = enricher(darkgrey_genes, TERM2GENE=term2gene)
head(summary(a))

# Save gene over representation analysis table
write.table(a, file = "LFMM_Darkgrey_Gene_over-representation.txt", sep = "\t", quote = F)

p7 <- barplot(a, main="Over Representation Analysis for LFMM Darkgrey Module") +
  ggtitle("Gene Over Representation Analysis for LFMM Darkgrey Module")

# Save the ggplot as a .tiff file
ggsave("LFMM_Darkgrey_gene_over_representation_barplot.tiff", plot = p7, device = "tiff", width = 7, height = 4)

p8 <- dotplot(a) +
  ggtitle("Gene Over Representation Analysis for LFMM Darkgrey Module")
ggsave("LFMM_Drakgrey_gene_over_representation_dotplot.tiff", plot = p8, device = "tiff", width = 7, height = 4)

#------------------------------------------------------------
# Enrichment analysis for genes LFMM darturqouise module
b = enricher(darkturqouise_genes, TERM2GENE=term2gene)
head(summary(b))

# Save gene over representation analysis table
write.table(b, file = "LFMM_Darkturqouise_Gene_over-representation.txt", sep = "\t", quote = F)

p9 <- barplot(b, main="Over Representation Analysis for LFMM Darkturqouise Module") +
  ggtitle("Gene Over Representation Analysis for LFMM Darkturqouise Module")

# Save the ggplot as a .tiff file
ggsave("LFMM_Darkturqouise_gene_over_representation_barplot.tiff", plot = p9, device = "tiff", width = 7, height = 4)

p10 <- dotplot(b) +
  ggtitle("Gene Over Representation Analysis for LFMM Darkturqouise Module")
ggsave("LFMM_Darkturqouise_gene_over_representation_dotplot.tiff", plot = p10, device = "tiff", width = 7, height = 4)

#------------------------------------------------------------
# Enrichment analysis for genes LFMM green module
c = enricher(green_genes, TERM2GENE=term2gene)
head(summary(c))

# Save gene over representation analysis table
write.table(c, file = "LFMM_Green_Gene_over-representation.txt", sep = "\t", quote = F)

p11 <- barplot(c, main="Over Representation Analysis for LFMM Green Module") +
  ggtitle("Gene Over Representation Analysis for LFMM Green Module")

# Save the ggplot as a .tiff file
ggsave("LFMM_Green_gene_over_representation_barplot.tiff", plot = p11, device = "tiff", width = 7, height = 4)

p12 <- dotplot(c) +
  ggtitle("Gene Over Representation Analysis for LFMM Green Module")
ggsave("LFMM_Green_gene_over_representation_dotplot.tiff", plot = p12, device = "tiff", width = 7, height = 4)


#------------------------------------------------------------
# Enrichment analysis for eGenes
d = enricher(eGene_genes, TERM2GENE=term2gene)
head(summary(d))

# Save gene over representation analysis table
write.table(d, file = "eGene_over-representation.txt", sep = "\t", quote = F)

p13 <- barplot(d, main="Over Representation Analysis for eGenes") +
  ggtitle("Gene Over Representation Analysis for eGenes")

# Save the ggplot as a .tiff file
ggsave("eGene_over_representation_barplot.tiff", plot = p13, device = "tiff", width = 7, height = 4)

p14 <- dotplot(d) +
  ggtitle("Gene Over Representation Analysis for eGenes")

ggsave("eGene_over_representation_dotplot.tiff", plot = p14, device = "tiff", width = 7, height = 4)


#------------------------------------------------------------
# We only performed gene overreprentation analysis for our selection outlier genes

# But here is how to go about gene set enrichment analysis if needed
#------------------------------------------------------------
#GSEA
# Load gene expression data
# df = read.table("./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/Deseq_result.txt", header=TRUE)

# load selection outliers with HAN412 gene IDs and corresponding ncbi gene ID
# outlier_genes <- read.table("./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/outlier_ncbi_blast.txt")
# names(outlier_genes) <- c("genes", "ncbi_id")

# subset expression dataset to contain only outlier genes
# df <- df %>%
# filter(genes %in% outlier_genes$genes)

# add ncbi gene id column to expression dataframe
# df_merge <- merge(df, outlier_genes, by = "genes", 
#                 all.x = TRUE) 


# we want the log2 fold change 
# original_gene_list <- df_merge$log2FoldChange

# name the vector
# names(original_gene_list) <- df_merge$ncbi_id

# omit any NA values 
# gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
# gene_list = sort(gene_list, decreasing = TRUE)

# convert transcript ID to ncbi loc ID
# gff_file <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/all_genes_exons.gff"

# Import the GFF3 file
# gff_data <- import(gff_file)

# Extract the ontology column (column 9) from the imported object
# gff_data2 <- data.frame(gff_data@elementMetadata)

#repair ID column in gff file 
# Extract the ID column from the element metadata
# ID <- gff_data2$ID

# Remove prefixes from the ID column
# ID <- gsub("mRNA-", "", ID)
# ID <- gsub("gene-", "", ID)
# ID <- gsub("cds-", "", ID)
# ID <- gsub("exon-", "", ID)
# ID <- gsub("rna-", "", ID)
# gff_data2$ID <- ID

# select IDs with corresponding loc names from gff file
# outlier_genes2 <- gff_data2 %>%
# filter(ID %in% names(gene_list)) %>%
# dplyr::select(., ID, gene) %>%
# distinct()

# gene_list2 <- do.call(rbind, lapply(gene_list, as.data.frame))
# gene_list3 <- tibble::rownames_to_column(gene_list2, "ID")

# gene_list4 <- merge(gene_list3, outlier_genes2, by = "ID", 
#                all.x = TRUE) 
# names(gene_list4) <- c("ID", "log2Change", "gene")

# final_gene_list <- gene_list4 %>%
# dplyr::select(., log2Change, gene) %>%
# na.omit() %>%
# distinct()

# Convert final genelist to named vector
# final_list2 <- final_gene_list %>%
# pull(log2Change, gene)

# sort the list in decreasing order (required for clusterProfiler)
# final_list2 = sort(final_list2, decreasing = TRUE)

# Get another term2gene for only genes in my expression dataset
# all_genes2 <- read.delim("./GCF_002127325.2_HanXRQr2.0-SUNRISE_gene_ontology.gaf", skip = 8, header = TRUE)

# subset all_genes2 to contain only genes in my expression dataset
# all_genes3 <- all_genes2 %>%
# filter(Symbol %in% names(final_list2)) %>%
# dplyr::select(., Symbol, GO_ID) %>%
# distinct()

# Get GO description for all genes2
# go_terms <- goIdToTerm(all_genes3$GO_ID)
# go_terms2 <- data.frame(go_terms)
# go_terms2$GO_ID <- names(go_terms)

# all_genes3$GO_terms <- go_terms2$go_terms

# term2gene2=all_genes3[, c("GO_terms", "Symbol")] %>%
# na.omit

# y <- GSEA(final_list2, TERM2GENE=term2gene2, pvalueCutoff = 1)
# head(summary(y))
# write.table(y, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Gene_set_enrichment.txt", sep = "\t", quote = F)

# p3 <- dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign) +
# ggtitle("GSEA")
# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/GSEA_differential_expression.tiff", plot = p3, device = "tiff", width = 15, height = 10)

# p4 <- dotplot(y) +
# ggtitle("GSEA")
# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/GSEA_differential_expression_enriched_ontology.tiff", plot = p4, device = "tiff")


# p5 <- ridgeplot(y) + 
# theme(axis.text.y=element_text(size=7),
#    axis.text.x=element_text(size=7)) +
# labs(x = "enrichment distribution")
# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/GSEA_enrichment_distribution.tiff", plot = p5, device = "tiff", width = 15, height = 10)

#---------------------------------------------------------------------------------

#KEGG
# search_kegg_organism('han', by='kegg_code')

# han <- search_kegg_organism('Helianthus annuus', by='scientific_name')
# dim(han)

# Get ncbi gene IDs for LFMM and PCAdapt outliers
# lfmm_gff <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/lfmm_exons.gff"

# Read the GFF file
# lfmm_gff <- readGFF(lfmm_gff)

# Extract gene IDs
# lfmm_genes <- lfmm_gff$gene %>% 
# str_replace("LOC","") %>%
# unique()

# pcadapt_gff <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/pcadapt_exons.gff"

# Read the GFF file
# pcadapt_gff <- readGFF(pcadapt_gff)

# Extract gene IDs
# pcadapt_genes <- pcadapt_gff$gene %>% 
# str_replace("LOC","") %>%
# unique()

# outlier_genes <- append(lfmm_genes, pcadapt_genes)

# KEGG over representation analysis for lfmm outlier genes
# kk <- enrichKEGG(gene         = lfmm_genes,
#               organism     = 'han',
#              pvalueCutoff = 0.05)
# head(kk) 
# write.table(kk, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/LFMM_KEGG_over-representation.txt", sep = "\t", quote = F)

# p6 <- barplot(kk) +
# ggtitle("KEGG Pathway Over Representation Analysis for LFMM Outlier Genes")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/LFMM_KEGG_pathway_OVA.tiff", plot = p6, device = "tiff", width = 7, height = 4)

# p7 <- dotplot(kk) +
# ggtitle("KEGG Pathway Over Representation Analysis for LFMM Outlier Genes")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/LFMM_KEGG_pathway_OVA_dotplot.tiff", plot = p7, device = "tiff", width = 7, height = 4)

# pathway <- data.frame(kk)

# Save KEGG ORA data
# write.table(pathway, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/LFMM_KEGG_over-representation.txt", sep = "\t", quote = F)

# browseKEGG(kk, 'han04261')

#------------------------------------------------------------------------------------
# KEGG Overrepresentation analysis for PCAdapt outlier genes
# kk2 <- enrichKEGG(gene         = pcadapt_genes,
#               organism     = 'han',
#              pvalueCutoff = 0.05)
# head(kk2) 
# write.table(kk2, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/PCAdapt_KEGG_over-representation.txt", sep = "\t", quote = F)

# p8 <- barplot(kk2) +
# ggtitle("KEGG Pathway Over Representation Analysis for PCAdapt Outlier Genes")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/PCAdapt_KEGG_pathway_OVA.tiff", plot = p8, device = "tiff", width = 7, height = 4)

# p9 <- dotplot(kk2) +
# ggtitle("KEGG Pathway Over Representation Analysis PCAdapt Outlier Genes")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/PCAdapt_KEGG_pathway_OVA_dotplot.tiff", plot = p9, device = "tiff", width = 7, height = 4)

# pathway <- data.frame(kk2)

# Save KEGG ORA data
# write.table(pathway, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/PCAdapt_KEGG_over-representation.txt", sep = "\t", quote = F)

# browseKEGG(kk2, 'han05168')

#------------------------------------------------------------------------------------




#KEGG pathway enrichment analysis
# change Loc IDs to NCBI IDs
# final_gene_list$gene <- gsub("^.{0,3}", "", final_gene_list$gene)

# Convert final genelist to named vector
# final_list3 <- final_gene_list %>%
# pull(log2Change, gene)

# sort the list in decreasing order (required for clusterProfiler)
# final_list3 = sort(final_list3, decreasing = TRUE)


# kk2 <- gseKEGG(geneList     = final_list3,
#             organism     = 'han',
#            minGSSize    = 120,
#           pvalueCutoff = 0.05,
#          verbose      = FALSE)
# head(kk2)
# write.table(kk2, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_enrichment.txt", sep = "\t", quote = F)

# p8 <- dotplot(kk2, showCategory=10, split=".sign") + 
# facet_grid(.~.sign) +
# ggtitle("KEGG Pathway GSEA")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_GSEA_differential_expression.tiff", plot = p8, device = "tiff")

# p9 <- dotplot(kk2) +
# ggtitle("KEGG Pathway GSEA")

# ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_GSEA_differential_expression.tiff_dot_plot", plot = p8, device = "tiff")

# Produce the native KEGG plot (PNG)
# dme <- pathview(gene.data=final_list3, pathway.id="han00053", species = 'han')

# Produce a different plot (PDF) (not displayed here)
# dme <- pathview(gene.data=final_list3, pathway.id="han00053", species = 'han', kegg.native = F)


# pathway for GSEA KEGG
# dme2 <- pathview(gene.data=final_list3, pathway.id="han01100", species = 'han')

# Produce a different plot (PDF) (not displayed here)
# dme2 <- pathview(gene.data=final_list3, pathway.id="han01100", species = 'han', kegg.native = F)

# browseKEGG(kk, 'han01100')
