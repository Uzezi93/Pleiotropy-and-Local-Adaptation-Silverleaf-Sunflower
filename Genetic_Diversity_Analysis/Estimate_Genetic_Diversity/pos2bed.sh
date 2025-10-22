#!/usr/bin/env bash

# A script to extract bed files from all test set positions 
set -euo pipefail

POS="eGenes_pos.txt"
GENE_BED="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/gene_coord.bed.gz"
OUT="egenes_coords.bed"

# load bedtools
module load bedtools2/2.31.1 || true

# Make a 0-based, half-open BED from 1-based single positions:
#   chrom  start(=pos-1)  end(=pos)
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' "$POS" \
  | LC_ALL=C sort -k1,1 -k2,2n > pos_tmp.bed

# Intersect: keep unique gene intervals from gene_coord that overlap any position
# gene_coord.bed.gz columns: chrom  start0  end1  gene  .  strand
bedtools intersect -a "$GENE_BED" -b pos_tmp.bed -wa -u \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -u > "$OUT"

rm -f pos_tmp.bed

echo "[OK] wrote $OUT"
