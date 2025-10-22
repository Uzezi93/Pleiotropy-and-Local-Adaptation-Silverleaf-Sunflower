#!/usr/bin/env bash
#SBATCH -J fc_array
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -p gpu
#SBATCH -G 8
#SBATCH -t 24:00:00
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=END

set -euo pipefail

# =================== Configuration ===================
GTF=${GTF:-genome/helianthus.gtf}
BAM_DIR=${BAM_DIR:-rg_bam}
OUT_DIR=${OUT_DIR:-counts_by_sample}
#STRAND=${STRAND:-2}              # 0=unstranded, 1=forward, 2=reverse
#MAPQ=${MAPQ:-10}                 # mapping-quality floor
MAX_PARALLEL=${MAX_PARALLEL:-20} # max concurrent array tasks
SUBMIT_COMBINE=${SUBMIT_COMBINE:-false}  # true => submit combine step automatically
MODE=${MODE:-work}               # internal: "work" or "combine"
# ================================================================================

# Activate conda env WITHOUT sourcing ~/.bashrc
if [ -f "$HOME/.conda/etc/profile.d/conda.sh" ]; then
  source "$HOME/.conda/etc/profile.d/conda.sh"
  conda activate myenv || true
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
  source "/opt/conda/etc/profile.d/conda.sh"
  conda activate myenv || true
else
  echo "[WARN] conda.sh not found; assuming featureCounts is on PATH"
fi

mkdir -p "$OUT_DIR"

# Build list of BAMs
BAM_LIST="$BAM_DIR/bams.list"
ls -1 "$BAM_DIR"/*.rg.md.bam | sort > "$BAM_LIST" || true
NBAMS=$(wc -l < "$BAM_LIST" 2>/dev/null || echo 0)

# ----------------------------- Launcher path ------------------------------------
if [[ "${MODE}" = "work" && -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ "$NBAMS" -lt 1 ]]; then
    echo "[ERROR] No BAMs found in $BAM_DIR (*.rg.md.bam)."
    exit 1
  fi
  SCRIPT=$(readlink -f "$0" 2>/dev/null || realpath "$0" 2>/dev/null || echo "$0")

  # Submit the worker array
  CHILD_JID=$(sbatch --parsable --array=1-"$NBAMS"%${MAX_PARALLEL} "$SCRIPT")
  echo "[INFO] Submitted array $CHILD_JID (1-$NBAMS%$MAX_PARALLEL)"

  # combine step after array finishes
  if [[ "$SUBMIT_COMBINE" == "true" ]]; then
    COMB_JID=$(sbatch --parsable --dependency=afterok:$CHILD_JID --export=ALL,MODE=combine "$SCRIPT")
    echo "[INFO] Submitted combine job $COMB_JID (after $CHILD_JID)"
  fi
  exit 0
fi

# ----------------------------- Combine path -------------------------------------
if [[ "${MODE}" = "combine" ]]; then
  echo "[INFO] Combining per-sample counts in $OUT_DIR -> counts_all.tsv"
  mapfile -t files < <(ls -1 "$OUT_DIR"/*.counts 2>/dev/null | sort)
  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "[ERROR] No *.counts files to combine in $OUT_DIR"; exit 1
  fi

  cut -f1 "${files[0]}" > "$OUT_DIR/.geneids.ref"
  for f in "${files[@]}"; do
    diff -q "$OUT_DIR/.geneids.ref" <(cut -f1 "$f") >/dev/null || {
      echo "[ERROR] Gene ID order mismatch among *.counts; aborting."
      exit 1
    }
  done

  paste "${files[@]}" \
  | awk 'BEGIN{OFS="\t"}
         NR==1 { printf $1; for (i=2;i<=NF;i+=2) printf "\t" $i; print ""; next }
               { printf $1; for (i=2;i<=NF;i+=2) printf "\t" $i; print "" }' \
  > counts_all.tsv

  rm -f "$OUT_DIR/.geneids.ref"
  echo "[INFO] Wrote counts_all.tsv"
  exit 0
fi

# ----------------------------- array  --------------------------------------
# Validate array index
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[ERROR] Worker mode requires SLURM_ARRAY_TASK_ID"; exit 1
fi
TASK_ID=${SLURM_ARRAY_TASK_ID}
if [[ "$TASK_ID" -lt 1 || "$TASK_ID" -gt "$NBAMS" ]]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$TASK_ID out of range (1..$NBAMS)"; exit 1
fi

BAM=$(sed -n "${TASK_ID}p" "$BAM_LIST")
base=$(basename "$BAM")
SAMPLE=${base%.rg.md.bam}

echo "[INFO] Task $TASK_ID/$NBAMS  Sample: $SAMPLE"
echo "[INFO] BAM: $BAM"
echo "[INFO] GTF: $GTF"

# Run featureCounts (paired-end)
featureCounts \
  -T "${SLURM_CPUS_PER_TASK:-8}" \
  -a "$GTF" \
  -o "${OUT_DIR}/${SAMPLE}.featureCounts.tsv" \
  -p -B -C \
  "$BAM"

# Make compact two-column file: Geneid \t <SAMPLE>
printf "Geneid\t%s\n" "$SAMPLE" > "${OUT_DIR}/${SAMPLE}.counts"
awk 'BEGIN{OFS="\t"} $1!="Geneid" && $1!~/^#/ {print $1, $NF}' \
  "${OUT_DIR}/${SAMPLE}.featureCounts.tsv" >> "${OUT_DIR}/${SAMPLE}.counts"

echo "[INFO] Wrote: ${OUT_DIR}/${SAMPLE}.featureCounts.tsv and ${OUT_DIR}/${SAMPLE}.counts"
