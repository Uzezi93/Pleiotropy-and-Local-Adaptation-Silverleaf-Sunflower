#!/usr/bin/env bash
#SBATCH -J divide_vcf
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 48:00:00
#SBATCH -p cpu
#SBATCH -o slurm-%A_%a.out
#SBATCH --array=1-6            # <-- one task per BED below

# Subset invariant vcf file into the different test sets

set -euo pipefail

module load bcftools/1.19
module load vcftools/0.1.16

# -------- inputs --------
VCF="../../../pixy_diploidized.vcf.gz"
OUTDIR="vcf_by_gene_sets"
mkdir -p "${OUTDIR}"

# List of BEDs to extract (order matches the --array range)
BEDS=(
  "eGene_genes.coords.bed"
  "eQTL_genes.coords.bed"
  "control_genes.coords.bed"
  "lfmm_genes.coords.bed"
  "pcadapt_genes.coords.bed"
  "shared_genes.coords.bed"
)

# -------- pick task's BED --------
i=${SLURM_ARRAY_TASK_ID}
BED=${BEDS[$((i-1))]}

if [[ ! -s "$BED" ]]; then
  echo "[ERROR] BED not found or empty: $BED" >&2
  exit 2
fi
if [[ ! -s "$VCF" ]]; then
  echo "[ERROR] VCF not found: $VCF" >&2
  exit 3
fi

# derive an output stem from the BED name
stem=$(basename "$BED")
stem=${stem%.bed}
stem=${stem%.bed.gz}
stem=${stem%.coords}          # e.g., eQTL_genes.coords.bed -> eQTL_genes

echo "[INFO] Task $i extracting with BED: $BED -> ${OUTDIR}/${stem}.vcf.gz"

# vcftools write to stdout, then bgzip + index
vcftools \
  --gzvcf "$VCF" \
  --bed "$BED" \
  --recode --recode-INFO-all --stdout \
| bgzip -c > "${OUTDIR}/${stem}.vcf.gz"

tabix -f -p vcf "${OUTDIR}/${stem}.vcf.gz"

echo "[DONE] ${OUTDIR}/${stem}.vcf.gz"
