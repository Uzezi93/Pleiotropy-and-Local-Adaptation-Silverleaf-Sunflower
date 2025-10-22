#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 10:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

# Script to extract positions of test sets from their respective vcf files, for sebequenctly getting their bed file positions

set -euo pipefail

shopt -s nullglob

# List your VCFs here (bgzipped + indexed)
vcfs=(
  eGenes.vcf.gz
  eQTL_genes.vcf.gz
  lfmm_genes.vcf.gz
  shared_genes.vcf.gz
  pcadapt_genes.vcf.gz
  control_genes.vcf.gz
)

# 1) Per-file outputs: <name>.sites.tsv (CHROM<TAB>POS), sorted & unique
for v in "${vcfs[@]}"; do
  out="${v%.vcf.gz}.sites.tsv"
  echo "[make] $out"
  bcftools query -f '%CHROM\t%POS\n' "$v" \
    | sort -k1,1 -k2,2n -u > "$out"
done

# 2) One combined file with a 'set' column
echo -e "set\tCHROM\tPOS" > all_sites.tsv
for v in "${vcfs[@]}"; do
  setname="${v%.vcf.gz}"
  bcftools query -f "${setname}\t%CHROM\t%POS\n" "$v" >> all_sites.tsv
done

echo "[done] wrote per-file *.sites.tsv and all_sites.tsv"
