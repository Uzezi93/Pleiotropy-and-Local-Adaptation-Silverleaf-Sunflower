#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

# Annotate vcf file with gene information to identify genes housing adaptive SNPs

module load bcftools/1.15

bcftools annotate \
  -a gene_coord.bed.gz \
  -c CHROM,FROM,TO,GENE \
  -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  filtered_snp.recode.vcf.gz > annotated.vcf
