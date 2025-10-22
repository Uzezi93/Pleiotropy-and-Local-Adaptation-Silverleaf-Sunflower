#!/bin/bash
#SBATCH -J vcf_merge
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -o slurm-%j.out
#SBATCH -G 4

# Recoding vcf file to 012 format for selection tests
plink --vcf filtered_snps.subset.vcf.gz --double-id --keep-allele-order --allow-extra-chr --make-bed --out argo_snps
plink --bfile argo_snps --recode A --allow-extra-chr --out argo_snps

# Output: argo_snps.raw  (0/1/2 counts for A1 allele per SNP; missing = NA)
