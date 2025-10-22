#!/bin/bash
#SBATCH -J pixy
#SBATCH -c 2
#SBATCH --mem=80G
#SBATCH -t 48:00:00
#SBATCH -p cpu
#SBATCH -o slurm-%j.out

# Load your conda environment

module load conda/latest

conda activate my_env

module load bcftools/1.19

# create a filtered VCF containing only invariant sites

vcftools --gzvcf variants2/chunks/all_filtered_snps.subset.vcf.gz \
 --max-maf 0 \
 --recode --stdout | bgzip -c > invariant.vcf.gz


 vcftools --gzvcf invariant.vcf.gz \
 --remove-indels \
 --recode --stdout | bgzip -c > invariant_snps.vcf.gz


# create a filtered VCF containing only variant sites
 vcftools --gzvcf variants2/chunks/all_filtered_snps.subset.vcf.gz \
 --mac 1 \
 --recode --stdout | bgzip -c > variant.vcf.gz

# index both vcfs using tabix
 bcftools index -t invariant_snp.vcf.gz
 tabix variant.vcf.gz
 tabix variants/chunk/filtered_snps.subset.vcf.gz # Used my filtered biallelic snp set instead

# combine the two VCFs using bcftools concat
# 1) Keep only rows with ALT='.'  (invariant ref sites)
 bcftools view -i 'ALT="."' invariant.vcf.gz -Oz -o invariant_refonly.vcf.gz
 tabix -f -p vcf invariant_refonly.vcf.gz
 bcftools index -t -f invariant_refonly.vcf.gz


 bcftools concat \
 --allow-overlaps \
  variants/chunks/filtered_snps.subset.vcf.gz invariant_refonly.vcf.gz \
  -O z > pixy_filtered.vcf.gz

 bcftools index -t pixy_filtered.vcf.gz


