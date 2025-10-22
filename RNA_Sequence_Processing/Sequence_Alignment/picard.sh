#!/usr/bin/env bash
#SBATCH -J picard.rg.md
#SBATCH -c 8                       # CPUs per task
#SBATCH --mem=40G                  # RAM per task
#SBATCH -p gpu                     # your site needs gpu partition
#SBATCH -G 4                       # (STAR/Picard won't use GPUs; quota requirement)
#SBATCH -t 06:00:00
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=END
#SBATCH --array=0-26


set -euo pipefail

IN_DIR="alignments_STAR"
OUT_DIR="rg_bam"
SAMPLES_FILE="samples.list"

mkdir -p "$OUT_DIR"

# Build a samples list (basenames) if it doesn't exist
if [[ ! -s "$SAMPLES_FILE" ]]; then
  ls -1 ${IN_DIR}/*.Aligned.sortedByCoord.out.bam \
    | sed 's#.*/##; s/\.Aligned\.sortedByCoord\.out\.bam$//' > "$SAMPLES_FILE"
fi

# Load samples into array and select this task's sample
mapfile -t SAMPLES < "$SAMPLES_FILE"
: "${SLURM_ARRAY_TASK_ID:=0}"
if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#SAMPLES[@]} )); then
  echo "Array index $SLURM_ARRAY_TASK_ID out of range 0..$((${#SAMPLES[@]}-1))" >&2
  exit 1
fi
sample="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

bam_in="${IN_DIR}/${sample}.Aligned.sortedByCoord.out.bam"
bam_rg="${OUT_DIR}/${sample}.rg.bam"
bam_md="${OUT_DIR}/${sample}.rg.md.bam"
dup_metrics="${OUT_DIR}/${sample}.dup_metrics.txt"

# -------- Locate a Picard JAR (prefer env, then module, then search) --------
PICARD_JAR="${PICARD_JAR:-}"
if [[ -z "${PICARD_JAR}" || ! -f "$PICARD_JAR" ]]; then
  # Try module (this site’s module exports $PICARD pointing to the jar)
  module load picard/3.1.1 2>/dev/null || true
  if [[ -n "${PICARD:-}" && -f "${PICARD}" ]]; then
    PICARD_JAR="${PICARD}"
  fi
fi
if [[ -z "${PICARD_JAR}" || ! -f "$PICARD_JAR" ]]; then
  # Last resort: search Spack tree
  PICARD_JAR="$(find /modules/spack/packages -type f -name 'picard*.jar' 2>/dev/null | head -n1 || true)"
fi
if [[ -z "${PICARD_JAR}" || ! -f "$PICARD_JAR" ]]; then
  echo "ERROR: Could not find picard*.jar. Set PICARD_JAR or fix the picard module." >&2
  exit 2
fi

# -------- Java (use system java; Picard 3.x is fine with >=17) --------
JAVA="java -Xmx30g"
# Per-node temp (keeps picard temp off NFS and speeds up)
TMPDIR="${TMPDIR:-/tmp/picard_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${sample}}"
mkdir -p "$TMPDIR"

echo "[$(date)] Sample: $sample"
echo "Java : $($JAVA -version 2>&1 | head -n1)"
echo "JAR  : $PICARD_JAR"
echo "Input: $bam_in"
echo "TMP  : $TMPDIR"

# Sanity check input
if [[ ! -s "$bam_in" ]]; then
  echo "ERROR: Missing input BAM: $bam_in" >&2
  exit 3
fi

# 1) Add/Replace read groups (idempotent)
if [[ ! -s "$bam_rg" ]]; then
  $JAVA -Djava.io.tmpdir="$TMPDIR" -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$bam_in" \
    O="$bam_rg" \
    RGID="$sample" \
    RGLB="lib1" \
    RGPL="ILLUMINA" \
    RGPU="${sample}.unit1" \
    RGSM="$sample" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT
else
  echo "RG BAM exists, skipping: $bam_rg"
fi

# 2) Mark duplicates (idempotent)
if [[ ! -s "$bam_md" ]]; then
  $JAVA -Djava.io.tmpdir="$TMPDIR" -jar "$PICARD_JAR" MarkDuplicates \
    I="$bam_rg" \
    O="$bam_md" \
    M="$dup_metrics" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT
else
  echo "MD BAM exists, skipping: $bam_md"
fi

echo "[$(date)] Done: $sample"
