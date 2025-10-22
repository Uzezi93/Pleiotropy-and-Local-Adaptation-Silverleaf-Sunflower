#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p cpu  # Partition
#SBATCH -G  0 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load blast-plus/2.14.1
# module load seqtk/1.3

# Get fasta file for all genes in our expression data
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta expressed_genes.txt >  final_enrichment2/expressed_genes.fasta

# Get fasta file for outlier  genes
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta outlier_genes.txt >  final_enrichment2/outlier_genes.fasta

# Get fasta file for outlier  shared genes
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta common_outliers_genes2.txt  >  final_enrichment2/shared_genes.fasta

# Get fasta file for genes in blue  module
# seqtk/seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta genes_in_blue.txt >  final_enrichment2/blue_genes.fasta

# Get fasta file for genes in green  module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta genes_in_green.txt > final_enrichment2/green_genes.fasta

# Get fasta file for genes in yellow module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta genes_in_yellow.txt >  final_enrichment2/yellow_genes.fasta

# Get fasta file for genes identified in pink module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta genes_in_pink.txt >  final_enrichment2/pink_genes.fasta

# Get fasta file for genes identified in red module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta genes_in_red.txt >  final_enrichment2/red_genes.fasta

# Get fasta file for eGenes
seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta eGenes.txt >  final_enrichment2/eGenes.fasta


# Create a BLAST database from the sunflower genome FASTA file
# makeblastdb -in NCBI/GCF_002127325.2_HanXRQr2.0-SUNRISE_cds_from_genomic.fa -dbtype nucl -out sunflower_db

# Blast the module fasta files  query against the sunflower genome database
# for i in final_enrichment2/*.fasta; do
  # blastn -query "$i" -db sunflower_db -outfmt "6 qseqid sacc" -out "final_enrichment2/blast_result/$(basename "$i").txt"
# done

