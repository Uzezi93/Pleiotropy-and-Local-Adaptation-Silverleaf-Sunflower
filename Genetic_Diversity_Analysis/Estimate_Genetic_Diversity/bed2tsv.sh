#!/usr/bin/env bash
set -euo pipefail

# Convert bed files to tsv before running angsd
# List of BED files
beds=(
  control_genes.coords.bed
  eGene_genes.coords.bed
  eQTL_genes.coords.bed
  lfmm_genes.coords.bed
  pcadapt_genes.coords.bed
  shared_genes.coords.bed
)

for bed in "${beds[@]}"; do
  [[ -s "$bed" ]] || { echo "[WARN] missing/empty: $bed"; continue; }
  base="${bed%.coords.bed}"

  echo "[INFO] Normalizing to TABs: $bed"
  # Re-write as tab-delimited (treat any whitespace as field sep, output tabs)
  # Keep only first 4 columns (chr, start, end, gene); strand is optional.
  awk 'BEGIN{OFS="\t"} NF>=4 {print $1,$2,$3,$4}' "$bed" > "${base}.tmp.bed"

  # Build ANGSD regions TSV: CHR:START-END<TAB>gene   (convert 0-based BED to 1-based)
  awk 'BEGIN{OFS="\t"} 
       NF>=4 {
         chr=$1; start=$2+1; end=$3; gene=$4;
         print chr ":" start "-" end, gene
       }' "${base}.tmp.bed" > "${base}.genes.regions.tsv"

  # Optional: regions-only list for -rf
  cut -f1 "${base}.genes.regions.tsv" > "${base}.genes.rf.txt"

  rm -f "${base}.tmp.bed"

  # Quick tab check
  awk -F'\t' 'NR==1{n=NF} NF!=n{bad=1} END{exit bad}' "${base}.genes.regions.tsv" \
    && echo "[OK] ${base}.genes.regions.tsv is tab-delimited" \
    || echo "[WARN] ${base}.genes.regions.tsv has inconsistent columns"

  echo "[DONE] ${base}: $(wc -l < "${base}.genes.regions.tsv") genes"
done
