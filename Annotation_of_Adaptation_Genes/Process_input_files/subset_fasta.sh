#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END


# Define the input FASTA file and the list of genes
fasta_file="Ha412HOv2.0-20181130.transcripts.fasta"
gene_list="LFMM_genes.txt"

# Loop through each gene in the list
while IFS= read -r gene; do
    # Filter headers based on the gene name and extract corresponding sequences
    grep -w -A9999999 "$gene" "$fasta_file" | awk -v gene="$gene" 'BEGIN {print ">"gene} {if (substr($0,1,1)==">") {next} else {print}}'
done < "$gene_list" > LFMM_genes.fasta
