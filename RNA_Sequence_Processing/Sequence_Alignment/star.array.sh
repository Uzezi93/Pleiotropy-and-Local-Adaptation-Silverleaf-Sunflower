#!/bin/bash
#SBATCH -c 8                          # CPUs per sample
#SBATCH --mem=40G                     # RAM per sample
#SBATCH -p gpu                        # your cluster's working partition
#SBATCH -G 4                          # STAR doesn't use GPUs, but gpu partition may require this
#SBATCH -t 06:00:00
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=END
#SBATCH --array=0-26                  # 27 samples -> indices 0..27  (edit if your count changes)

set -euo pipefail

# --- Conda activation (so STAR is found) ---
eval "$(conda shell.bash hook)"
conda activate myenv

which STAR
STAR --version

GENOME_DIR="index"                    # STAR genome index
IN_DIR="trimmed_fastq"
OUT_DIR="alignments_STAR"
mkdir -p "$OUT_DIR"

# Build sample list (handles .fq and .fq.gz)
mapfile -t SAMPLES < <(
  ls "$IN_DIR"/*_1.trim.fq "$IN_DIR"/*_1.trim.fq.gz 2>/dev/null \
  | sed 's#.*/##' | sed -E 's/_1\.trim\.fq(\.gz)?$//' | sort -u
)

# Bounds check
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "No input FASTQs found in $IN_DIR matching *_1.trim.fq[.gz]" >&2
  exit 1
fi
if [[ ${SLURM_ARRAY_TASK_ID:-0} -ge ${#SAMPLES[@]} ]]; then
  echo "ARRAY index ${SLURM_ARRAY_TASK_ID} >= number of samples ${#SAMPLES[@]}" >&2
  exit 1
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Pick files and decompressor
if [[ -f "$IN_DIR/${SAMPLE}_1.trim.fq.gz" ]]; then
  R1="$IN_DIR/${SAMPLE}_1.trim.fq.gz"
  R2="$IN_DIR/${SAMPLE}_2.trim.fq.gz"
  READCMD=(--readFilesCommand zcat)
else
  R1="$IN_DIR/${SAMPLE}_1.trim.fq"
  R2="$IN_DIR/${SAMPLE}_2.trim.fq"
  READCMD=()
fi

echo "[$(date)] STAR sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

STAR \
  --runThreadN "$THREADS" \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$R1" "$R2" \
  "${READCMD[@]}" \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outFileNamePrefix "${OUT_DIR}/${SAMPLE}."
