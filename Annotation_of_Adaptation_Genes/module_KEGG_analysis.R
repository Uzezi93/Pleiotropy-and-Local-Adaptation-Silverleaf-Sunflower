setwd("/project/pi_brook_moyers_umb_edu/Uzezi_argo/enrichment_analysis/")
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

# Get ncbi gene IDs for LFMM and PCAdapt outliers
steelblue_gff <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/modules/gff_files/steelblue_edge.gff"

# Read the GFF file
steelblue_gff <- readGFF(steelblue_gff)

# Extract gene IDs
steelblue_gff_genes <- steelblue_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = steelblue_gff_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

q1 <- barplot(kk) +
  ggtitle("Steelblue Pathway KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Steelblue_Pathway_KEGG_pathway_OVA.tiff", plot = q1, device = "tiff")

q2 <- dotplot(kk) +
  ggtitle("Ivory Module KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Steelblue_Pathway_KEGG_pathway_OVA_dotplot.tiff", plot = q2, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_pathway_ORA_steelblue_attribute.txt", sep = "\t", quote = F)

#----------------------------------------------------------------------------
# Get ncbi gene IDs for LFMM and PCAdapt outliers
grey_gff <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/modules/gff_files/grey_attribute.gff"

# Read the GFF file
grey_gff <- readGFF(grey_gff)

# Extract gene IDs
grey_genes <- grey_gff$gene %>% 
  str_replace("LOC","") 

# KEGG over representation analysis
kk <- enrichKEGG(gene         = grey_genes,
                 organism     = 'han',
                 pvalueCutoff = 1)
head(kk) 

q1 <- barplot(kk) +
  ggtitle("Yellow Pathway KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Yellow_Pathway_KEGG_pathway_OVA.tiff", plot = q1, device = "tiff")

q2 <- dotplot(kk) +
  ggtitle("yellow Module KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Yellow_Pathway_KEGG_pathway_OVA_dotplot.tiff", plot = q2, device = "tiff")

pathway <- data.frame(kk)

write.table(pathway, file = "../Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_pathway_ORA_yellow_attribute.txt", sep = "\t", quote = F)





#MKEGG pathway enrichment analysis
# change Loc IDs to NCBI IDs
# Load gene expression data
df = read.table("./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/Deseq_result.txt", header=TRUE)

# load selection outliers with HAN412 gene IDs and corresponding ncbi gene ID
ivory_genes <- read.table("./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/ivory_blast_results.txt")
names(ivory_genes) <- c("genes", "ncbi_id")

# subset expression dataset to contain only outlier genes
df <- df %>%
  filter(genes %in% ivory_genes$genes)

# add ncbi gene id column to expression dataframe
df_merge <- merge(df, ivory_genes, by = "genes", 
                  all.x = TRUE) 


# we want the log2 fold change 
original_gene_list <- df_merge$log2FoldChange

# name the vector
names(original_gene_list) <- df_merge$ncbi_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# convert transcript ID to ncbi loc ID
gff_file <- "./GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.gff"

# Import the GFF3 file
gff_data <- import(gff_file)

# Extract the ontology column (column 9) from the imported object
gff_data2 <- data.frame(gff_data@elementMetadata)

#repair ID column in gff file 
# Extract the ID column from the element metadata
ID <- gff_data2$ID

# Remove prefixes from the ID column
ID <- gsub("mRNA-", "", ID)
ID <- gsub("gene-", "", ID)
ID <- gsub("cds-", "", ID)
ID <- gsub("exon-", "", ID)
ID <- gsub("rna-", "", ID)
gff_data2$ID <- ID

# get corresponding loc names from gff file
ivory_genes2 <- gff_data2 %>%
  filter(ID %in% names(gene_list)) %>%
  dplyr::select(., ID, gene) %>%
  distinct()

gene_list2 <- do.call(rbind, lapply(gene_list, as.data.frame))
gene_list3 <- tibble::rownames_to_column(gene_list2, "ID")

gene_list4 <- merge(gene_list3, ivory_genes2, by = "ID", 
                    all.x = TRUE) 
names(gene_list4) <- c("ID", "log2Change", "gene")

final_gene_list <- gene_list4 %>%
  dplyr::select(., log2Change, gene) %>%
  na.omit() %>%
  distinct()

# Convert final genelist to named vector
final_list2 <- final_gene_list %>%
  pull(log2Change, gene)

# sort the list in decreasing order (required for clusterProfiler)
final_list2 = sort(final_list2, decreasing = TRUE)


final_gene_list$gene <- gsub("^.{0,3}", "", final_gene_list$gene)

# Convert final genelist to named vector
final_list3 <- final_gene_list %>%
  pull(log2Change, gene)

# sort the list in decreasing order (required for clusterProfiler)
final_list3 = sort(final_list3, decreasing = TRUE)


Mkk2 <- gseMKEGG(geneList     = final_list3,
               organism     = 'han',
               pvalueCutoff = 1)
head(Mkk2)

p8 <- dotplot(kk2, showCategory=10, split=".sign") + 
  facet_grid(.~.sign) +
  ggtitle("KEGG Pathway GSEA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_GSEA_differential_expression.tiff", plot = p8, device = "tiff")

p9 <- dotplot(kk2) +
  ggtitle("KEGG Pathway GSEA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/KEGG_GSEA_differential_expression.tiff_dot_plot", plot = p8, device = "tiff")

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=final_list3, pathway.id="han00053", species = 'han')

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=final_list3, pathway.id="han00053", species = 'han', kegg.native = F)


# pathway for GSEA KEGG
dme2 <- pathview(gene.data=final_list3, pathway.id="han01100", species = 'han')

# Produce a different plot (PDF) (not displayed here)
dme2 <- pathview(gene.data=final_list3, pathway.id="han01100", species = 'han', kegg.native = F)

browseKEGG(kk, 'han01100')


#-------------------------------------------------------------------------------
# For Steelblue
# Get ncbi gene IDs for LFMM and PCAdapt outliers
steelblue_gff <- "./20240301_HanXRQr2.0-20180814.vs.Ha412HOv2.0-20181130/steelblue_exons.gff"

# Read the GFF file
steelblue_gff <- readGFF(steelblue_gff)

# Extract gene IDs
steelblue_genes <- steelblue_gff$gene %>% 
  str_replace("LOC","") %>%
  unique()

# KEGG over representation analysis
kk <- enrichMKEGG(gene         = steelblue_genes,
                  organism     = 'han',
                  pvalueCutoff = 0.05)
head(kk) 

q3 <- barplot(kk) +
  ggtitle("Steelblue Module KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Steelblue_Module_KEGG_pathway_OVA.tiff", plot = q3, device = "tiff")

q4 <- dotplot(kk) +
  ggtitle("Steelblue Module KEGG Pathway ORA")

ggsave("/project/pi_brook_moyers_umb_edu/Uzezi_argo/Helianthus_argophyllus_project/argo_plots/PLOS_genetics_figures/Steelblue_Module_KEGG_pathway_OVA_dotplot.tiff", plot = q4, device = "tiff")

pathway <- data.frame(kk)


