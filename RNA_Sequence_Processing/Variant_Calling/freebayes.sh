#!/usr/bin/env bash
#SBATCH -J fb_array
#SBATCH -c 8                         # CPUs per array task (used by GNU parallel)
#SBATCH --mem=32G                    # RAM per task
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -t 24:00:00
#SBATCH -o slurm-%A_%a.freebayes.out
#SBATCH --array=0-99

set -euo pipefail

# --- modules (adjust to your cluster names/versions) ---
# module load freebayes/1.3.5
# module load parallel
# module load vcflib          # vcffirstheader, vcfstreamsort, vcfuniq
# module load htslib          # bgzip, tabix
# module load samtools

# --- paths ---
REF="genome/Ha412HOv2.0-20181130.fasta"
BAM_DIR="rg_bam"
BAM_LIST="${BAM_DIR}/bams.list"
REGIONS="genome/Ha412HOv2.fa.100kb.regions"

OUT_DIR="variants"
CHUNK_DIR="${OUT_DIR}/chunks"
mkdir -p "$OUT_DIR" "$CHUNK_DIR"

# --- threads per task ---
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- checks ---
for x in freebayes parallel vcffirstheader vcfstreamsort vcfuniq bgzip tabix samtools; do
  command -v "$x" >/dev/null || { echo "ERROR: missing tool '$x'"; exit 1; }
done

[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
[[ -s "$BAM_LIST" ]] || ls -1 ${BAM_DIR}/*.rg.md.bam > "$BAM_LIST"
[[ -s "$REGIONS" ]] || fasta_generate_regions.py "${REF}.fai" 100000 > "$REGIONS"

# --- compute array ---
MIN=${SLURM_ARRAY_TASK_MIN:-0}
MAX=${SLURM_ARRAY_TASK_MAX:-0}
STEP=${SLURM_ARRAY_TASK_STEP:-1}
IDX=${SLURM_ARRAY_TASK_ID:-0}
COUNT=$(( (MAX - MIN)/STEP + 1 ))
SLOT=$(( (IDX - MIN)/STEP ))

# --- pick a task ---
TASK_REGIONS="${CHUNK_DIR}/regions.${IDX}.list"
awk -v mod="$COUNT" -v rem="$SLOT" '((NR-1)%mod)==rem' "$REGIONS" > "$TASK_REGIONS"

# --- build FreeBayes base command ---
fb_cmd=( freebayes
  -f "$REF"
  -L "$BAM_LIST"
  --genotype-qualities
  --min-base-quality 10
  --min-mapping-quality 1
)

OUT_CHUNK="${CHUNK_DIR}/freebayes.chunk_${IDX}.vcf.gz"

echo "[$(date)] Task ${SLURM_ARRAY_JOB_ID}_${IDX}: regions=$(wc -l < "$TASK_REGIONS") threads=$THREADS"
cat "$TASK_REGIONS" \
  | parallel -k -j "$THREADS" "${fb_cmd[@]}" --region {} \
  | vcffirstheader \
  | vcfstreamsort -w 1000 \
  | vcfuniq \
  | bgzip -c > "$OUT_CHUNK"

tabix -p vcf "$OUT_CHUNK"
echo "[$(date)] Task ${SLURM_ARRAY_JOB_ID}_${IDX} done → $OUT_CHUNK"
