#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load vcftools/0.1.14
module load bcftools/1.15

# create new vcf subsets for eQTLs, eGenes, and Control to get genes for connectivity plot

vcftools --gzvcf annotated.vcf --positions eQTL_position2.txt --recode --out eQTL2_annotated

# Turn eGenes positions text file to tab delimited text file
awk -v OFS='\t' '{$1=$1}1' eGenes_position3.txt > eGenes_position4.txt

vcftools --gzvcf annotated.vcf --positions eGenes_position4.txt --recode --out eGenes2_annotated

vcftools --gzvcf annotated.vcf --exclude-positions combined_positions2.txt --recode --recode-INFO-all --out control2_annotated

bcftools query -f '%ID\n' control2_annotated.recode.vcf > control_genes.txt

# Get unique gene names for control genes
awk -F'[;]' '{sub(/^Name=/, ""); print $1}' control_genes.txt | uniq > control_genes2.txt

bcftools query -f '%ID\n' eGenes2_annotated.recode.vcf > eGenes_genes.txt

# Get unique gene names for eGenes
awk -F'[;]' '{sub(/^Name=/, ""); print $1}' eGenes_genes.txt | uniq > eGenes_genes2.txt

bcftools query -f '%ID\n' eQTL2_annotated.recode.vcf > eQTL_genes.txt

# Get unique eQTL gene names
awk -F'[;]' '{sub(/^Name=/, ""); print $1}' eQTL_genes.txt | uniq > eQTL_genes2.txt

# Get lfmm and pcadapt annotated vcf file
vcftools --gzvcf annotated.vcf --positions lfmm3_pos.txt --recode --out lfmm_annotated

vcftools --gzvcf annotated.vcf --positions pcadapt_pos.txt --recode --out pcadapt_annotated

