#!/bin/bash
#SBATCH -J vcf_merge
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -o slurm-%j.out
#SBATCH -G 4


set -euo pipefail

# Modules
module load bcftools/1.19 || true

# I also used this script for generating and filtering biallelic variants
# INPUT
# IN_VCF="freebayes.merged.norm.vcf.gz"
IN_VCF="../../invariant_refonly.vcf.gz"

# Outputs
# RAW_SNPS="invariant_snp.vcf.gz"
FILT_SNPS="invariant_filtered_snps.vcf.gz"

# 1) Keep **biallelic SNPs only** (strict SNP selection)
#    -m2 -M2 : only biallelic
#    -v snps : only SNPs (excludes MNPs/complex)
# bcftools view -m2 -M2 -v snps -Oz -o "$RAW_SNPS" "$IN_VCF"
# bcftools index -t "$RAW_SNPS"

# 2) Apply site/genotype filters with vcftools
# NOTE: --hwe 0.8 is extremely stringent; most workflows use something like 1e-6.
# Adjust to your study design. Here we set 1e-6 (change if you intend a different threshold).
vcftools \
  --gzvcf "$IN_VCF" \
  --minQ 30 \
  --minGQ 20 \
  --max-missing 0.8 \
  --min-meanDP 5 \
  --max-meanDP 40 \
  --recode --recode-INFO-all --stdout \
  | bgzip -c > "$FILT_SNPS"

# Index final
bcftools index -t "$FILT_SNPS"

echo "Done:"
# echo "  Biallelic SNPs : $RAW_SNPS (+ .tbi)"
echo "  Filtered SNPs  : $FILT_SNPS (+ .tbi)"
