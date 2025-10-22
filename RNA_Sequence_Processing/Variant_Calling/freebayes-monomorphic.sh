#!/usr/bin/env bash
#SBATCH -J fb_array
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH -t 48:00:00
#SBATCH -o slurm-%A_%a.freebayes.out
#SBATCH --array=0-99

set -euo pipefail

# --- modules ---
# module load freebayes/1.3.5
# module load parallel
# module load vcflib          # vcffirstheader, vcfstreamsort, vcfuniq
# module load htslib          # bgzip, tabix
# module load bcftools
# module load samtools

# --- paths ---
REF="genome/Ha412HOv2.0-20181130.fasta"
BAM_DIR="rg_bam"
BAM_LIST="${BAM_DIR}/bams.list"
REGIONS="genome/Ha412HOv2.fa.100kb.regions"   # master 100 kb windows

OUT_DIR="variants2"
CHUNK_DIR="${OUT_DIR}/chunks"
mkdir -p "$OUT_DIR" "$CHUNK_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- tool checks ---
for x in freebayes parallel vcffirstheader vcfstreamsort vcfuniq bgzip tabix bcftools samtools; do
  command -v "$x" >/dev/null || { echo "ERROR: missing tool '$x'"; exit 1; }
done

# --- ensure FASTA index & BAM list ---
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
[[ -s "$BAM_LIST" ]] || ls -1 ${BAM_DIR}/*.rg.md.bam > "$BAM_LIST"

# --- rebuild regions if missing (you deleted them) ---
if [[ ! -s "$REGIONS" ]]; then
  echo "[$(date)] Rebuilding regions @100kb → $REGIONS"
  fasta_generate_regions.py "${REF}.fai" 100000 > "$REGIONS"
fi

# --- array math ---
MIN=${SLURM_ARRAY_TASK_MIN:-0}
MAX=${SLURM_ARRAY_TASK_MAX:-0}
STEP=${SLURM_ARRAY_TASK_STEP:-1}
IDX=${SLURM_ARRAY_TASK_ID:-0}
COUNT=$(( (MAX - MIN)/STEP + 1 ))
SLOT=$(( (IDX - MIN)/STEP ))

# --- per-task region list (round-robin split) ---
TASK_REGIONS="${CHUNK_DIR}/regions.${IDX}.list"
awk -v mod="$COUNT" -v rem="$SLOT" '((NR-1)%mod)==rem' "$REGIONS" > "$TASK_REGIONS"

# --- FreeBayes command (emit variants + monomorphic) ---
fb_cmd=( freebayes
  -f "$REF"
  -L "$BAM_LIST"
  --genotype-qualities
  --min-base-quality 10
  --min-mapping-quality 1
  --min-coverage 1
  --use-reference-allele
  --report-monomorphic
)

OUT_ALL="${CHUNK_DIR}/freebayes.all_sites.chunk_${IDX}.vcf.gz"
OUT_MONO="${CHUNK_DIR}/freebayes.monomorphic_only.chunk_${IDX}.vcf.gz"
OUT_VAR="${CHUNK_DIR}/freebayes.variants_only.chunk_${IDX}.vcf.gz"

echo "[$(date)] Task ${SLURM_ARRAY_JOB_ID}_${IDX}: regions=$(wc -l < "$TASK_REGIONS") threads=$THREADS"
export PARALLEL="--halt soon,fail=1"

# 1) run per-region in parallel; consolidate; sort; uniq; compress
cat "$TASK_REGIONS" \
  | parallel -k -j "$THREADS" "${fb_cmd[@]}" --region {} \
  | vcffirstheader \
  | vcfstreamsort -w 1000 \
  | vcfuniq \
  | bgzip -c > "$OUT_ALL"

tabix -p vcf "$OUT_ALL"

# 2) split by ALT (monomorphic vs variant)
bcftools view -i 'ALT="."'  -O z -o "$OUT_MONO" "$OUT_ALL" && tabix -p vcf "$OUT_MONO"
bcftools view -i 'ALT!="."' -O z -o "$OUT_VAR"  "$OUT_ALL" && tabix -p vcf "$OUT_VAR"

# 3) quick counts
echo "records (all)  : $(bcftools index -n "$OUT_ALL")"
echo "records (mono) : $(bcftools index -n "$OUT_MONO")"
echo "records (var)  : $(bcftools index -n "$OUT_VAR")"

echo "[$(date)] Task ${SLURM_ARRAY_JOB_ID}_${IDX} done"
