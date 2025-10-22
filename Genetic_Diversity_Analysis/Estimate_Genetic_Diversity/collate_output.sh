#!/usr/bin/env bash
#SBATCH -J collate_angsd_pestpg
#SBATCH -c 24
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH -p cpu
#SBATCH -o slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --export=NONE

# Script to save all angsd output per gene into one csv file
set -euo pipefail
module purge || true              # clean env
export LC_ALL=C

# ---------- CONFIG ----------
ROOT="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets/angsd_per_gene_out"
OUT="${ROOT}/angsd_pestPG_collated_by_gene_pop.tsv"

# ---------- HEADER ----------
# Columns extracted from pestPG:
#   meta  Chr  WinCenter  tW  tP  tF  tH  tL  Tajima  fuf  fud  fayh  zeng  nSites
# We’ll keep: module  pop  gene  chr  win_center  tW  tP  Tajima  nSites  file
echo -e "module\tpop\tgene\tchromosome\twin_center\ttW\ttP\tTajima\tnSites\tfile" > "${OUT}"

shopt -s nullglob

# Loop genes.pop directories
for d in "${ROOT}"/*.*; do
  [[ -d "$d" ]] || continue
  base="$(basename "$d")"
  # module = everything before last dot; pop = after last dot
  module="${base%.*}"
  pop="${base##*.}"

  # Walk each gene's pestPG
  for f in "$d"/*.theta.thetas.idx.pestPG; do
    [[ -s "$f" ]] || continue

    # Gene tag: strip directory & suffix. My filenames look like "gene.<ID>.theta.thetas.idx.pestPG"
    fname="$(basename "$f")"
    gene="${fname%.theta.thetas.idx.pestPG}"
    # remove common "gene." prefix 
    gene="${gene#gene.}"

    # Parse pestPG: skip header lines starting with "#" and keep the lines starting with "("
    awk -v MOD="$module" -v POP="$pop" -v GENE="$gene" -v FILE="$f" -v OFS="\t" '
      BEGIN{ n=0 }
      $0 ~ /^#/ { next }                                # skip header comment line
      $1 ~ /^\(/ {
        chr=$2; win=$3; tW=$4; tP=$5; taj=$9; nSites=$14;
        # Coerce to numbers where possible (some files have trailing dots)
        gsub(/\.$/,"", nSites)
        print MOD, POP, GENE, chr, win, tW, tP, taj, nSites, FILE;
        n++
      }
      END{
        if(n==0){
          # No data rows; print an NA row to keep track
          # print MOD, POP, GENE, "NA", "NA", "NA", "NA", "NA", "0", FILE > "/dev/stderr"
        }
      }
    ' "$f" >> "${OUT}"
  done
done

echo "[OK] wrote: ${OUT}"
