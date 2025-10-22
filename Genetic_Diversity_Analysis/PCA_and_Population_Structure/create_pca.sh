#!/bin/bash
#SBATCH -J pca
#SBATCH -c 4
#SBATCH --mem=80G
#SBATCH -t 48:00:00
#SBATCH -p cpu
#SBATCH -o slurm-%j.out

# 1. Prune sites
plink2 --vcf variants/chunks/filtered_snps.subset.vcf \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --indep-pairwise 50 10 0.1 \
  --out argo

# 2. Create binary files (BED/BIM/FAM) for the full dataset
plink2 --vcf variants/chunks/filtered_snps.subset.vcf \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out argo

# 3. Calculate allele frequencies
plink2 --bfile argo \
  --allow-extra-chr \
  --freq \
  --out argo

# 4. Perform PCA using the frequency file
plink2 --bfile argo \
  --allow-extra-chr \
  --read-freq argo.afreq \
  --pca \
  --out argo
