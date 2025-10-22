#!/bin/bash
#SBATCH -J admixture_numeric_array
#SBATCH -c 10
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -G 8
#SBATCH -o slurm-%A_%a.out
#SBATCH --array=1-7                 # K = 1..7

set -euo pipefail

# ---- paths ----
SRC_PREFIX="../variants/chunks/argo_snps"     # my existing PLINK trio
NUM_PREFIX="../variants/chunks/argo_snps.admix"   # numeric copy to use for ADMIXTURE
OUTDIR="ADMIXTURE_REAL"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

# ---- required modules ----
module load admixture/1.3.0 2>/dev/null || module load admixture 2>/dev/null || true

# ---- sanity checks ----
for ext in bed bim fam; do
  [[ -s "${SRC_PREFIX}.${ext}" ]] || { echo "Missing ${SRC_PREFIX}.${ext}"; exit 2; }
done

# ---- Create a numeric BIM copy (maps ANY code, incl. 0/chr*/scaffolds, to 1..N) ----
echo "Creating numeric PLINK copy: ${NUM_PREFIX}.{bed,bim,fam}"
cp -f "${SRC_PREFIX}.bed" "${NUM_PREFIX}.bed"
cp -f "${SRC_PREFIX}.fam" "${NUM_PREFIX}.fam"

# map original chrom codes -> 1..N and save mapping
: > "${OUTDIR}/contig_map.tsv"
awk 'BEGIN{OFS="\t"}
{
  if(!($1 in m)){ m[$1]=++n; print m[$1], $1 >> "'"${OUTDIR}"'/contig_map.tsv" }
  $1=m[$1]; print
}' "${SRC_PREFIX}.bim" > "${NUM_PREFIX}.bim"

# quick check
echo "Original unique chrom codes: $(cut -f1 ${SRC_PREFIX}.bim | sort -u | wc -l)"
echo "Numeric unique chrom codes : $(cut -f1 ${NUM_PREFIX}.bim | sort -u | wc -l)"
echo "First few numeric BIM lines:"
head -3 "${NUM_PREFIX}.bim" || true

# ---- run ADMIXTURE on the NUMERIC copy ----
K=${SLURM_ARRAY_TASK_ID}
THREADS=${SLURM_CPUS_PER_TASK:-10}

echo "[$(date)] Running ADMIXTURE K=${K} on ${NUM_PREFIX}.bed with ${THREADS} threads"
admixture --cv=10 -j${THREADS} "${NUM_PREFIX}.bed" ${K} | tee "${LOGDIR}/K${K}.log"

# move outputs to a tidy location
mv -f "${NUM_PREFIX}.${K}.Q" "${OUTDIR}/K${K}.Q"
mv -f "${NUM_PREFIX}.${K}.P" "${OUTDIR}/K${K}.P"
awk '{print $2}' "${NUM_PREFIX}.fam" > "${OUTDIR}/samples.order"

# make a simple TSV to plot later
{
  printf "sample"
  for i in $(seq 1 ${K}); do printf "\tQ%d", i; done
  printf "\n"
  paste "${OUTDIR}/samples.order" "${OUTDIR}/K${K}.Q"
} > "${OUTDIR}/K${K}.Q.table.tsv"

# record CV error
grep -h "CV error" "${LOGDIR}/K${K}.log" \
  | sed -E 's/.*K=([0-9]+)\):[[:space:]]*([0-9.]+).*/\1\t\2/' \
  >> "${OUTDIR}/cv_summary.tsv"

echo "Done K=${K}. See ${OUTDIR}/K${K}.Q, ${OUTDIR}/K${K}.P, ${LOGDIR}/K${K}.log"
