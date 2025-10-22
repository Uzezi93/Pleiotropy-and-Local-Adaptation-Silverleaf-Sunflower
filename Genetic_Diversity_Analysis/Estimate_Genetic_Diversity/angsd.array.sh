#!/usr/bin/env bash
#SBATCH -J angsd_per_gene_multi
#SBATCH -c 20
#SBATCH --mem=16G
#SBATCH -t 48:00:00
#SBATCH -p gpu
#SBATCH -G 20
#SBATCH -o slurm-%A_%a.out
#SBATCH --array=1-1000%200

# Running ansgd per gene to compute Fay's H statistics
set -euo pipefail

# ==================== batching config ====================
# How many MANIFEST rows each array task should process
ITEMS_PER_TASK=40

# ==================== modules ====================
module load angsd/0.935
module load samtools/1.19.2 || true

# ==================== paths ====================
REF="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/genome/Ha412HOv2.0-20181130.fasta"
BAM_DIR="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/rg_bam"
WKDIR="$(pwd)"

# Population sample lists (plain text, one ID per line)
COAST_SAMPLES="${WKDIR}/coast.samples"
NORTH_SAMPLES="${WKDIR}/north.samples"

# Output layout
OUTROOT="${WKDIR}/angsd_per_gene_out"
BAMLIST_DIR="${OUTROOT}/bamlists"
LOGDIR="${OUTROOT}/logs"
MANIFEST="${OUTROOT}/task_manifest.tsv"   # columns: POP \t MODULE \t REGION \t GENE

mkdir -p "${OUTROOT}" "${BAMLIST_DIR}" "${LOGDIR}"

# ==================== helper: build bamlist ====================
build_bamlist () {
  local samples_file="$1"
  local bamlist_out="$2"
  : > "${bamlist_out}"

  while read -r id || [[ -n "${id}" ]]; do
    [[ -z "${id}" ]] && continue
    local md="${BAM_DIR}/${id}.rg.md.bam"
    local rg="${BAM_DIR}/${id}.rg.bam"
    local bam=""
    if   [[ -f "$md" ]]; then bam="$md"
    elif [[ -f "$rg" ]]; then bam="$rg"
    else
      echo "[WARN] No BAM for sample: ${id}" >&2; continue
    fi
    echo "$bam" >> "$bamlist_out"
    # ensure index
    [[ -f "${bam}.bai" ]] || { echo "[INFO] Indexing ${bam}"; samtools index -@ 2 "$bam" || true; }
  done < <(grep -v '^[[:space:]]*$' "${samples_file}")

  echo "[OK] Wrote $(wc -l < "${bamlist_out}") entries -> ${bamlist_out}"
}

# ==================== build BAM lists (once) ====================
[[ -s "${COAST_SAMPLES}" ]] || { echo "[ERR] Missing coast.samples" >&2; exit 2; }
[[ -s "${NORTH_SAMPLES}" ]] || { echo "[ERR] Missing north.samples" >&2; exit 3; }

COAST_BAMLIST="${BAMLIST_DIR}/coast.bamlist"
NORTH_BAMLIST="${BAMLIST_DIR}/north.bamlist"
[[ -s "${COAST_BAMLIST}" ]] || build_bamlist "${COAST_SAMPLES}" "${COAST_BAMLIST}"
[[ -s "${NORTH_BAMLIST}" ]] || build_bamlist "${NORTH_SAMPLES}" "${NORTH_BAMLIST}"

# ==================== discover & validate region files ====================
# Expect files: *.genes.regions.tsv (2 cols: REGION<TAB>GENE), REGION like CHR:START-END
mapfile -t REG_FILES < <(ls -1 "${WKDIR}"/*.genes.regions.tsv 2>/dev/null || true)
if (( ${#REG_FILES[@]} == 0 )); then
  echo "[ERR] no valid *.genes.regions.tsv files found in ${WKDIR}" >&2
  exit 6
fi

VALID_FILES=()
for f in "${REG_FILES[@]}"; do
  [[ -s "$f" ]] || { echo "[WARN] skipping empty: $f"; continue; }
  sed -i 's/\r$//' "$f"
  if awk -F'\t' 'NF && NF<2{bad=1} END{exit bad}' "$f"; then
    VALID_FILES+=("$f")
  else
    echo "[WARN] File not tab-delimited (2+ cols), skipping: $f"
  fi
done
[[ "${#VALID_FILES[@]}" -gt 0 ]] || { echo "[ERR] all *.genes.regions.tsv files are invalid" >&2; exit 6; }

# ==================== build manifest (once) ====================
# MANIFEST columns: POP \t MODULE \t REGION \t GENE
if [[ ! -s "${MANIFEST}" ]]; then
  : > "${MANIFEST}"
  for f in "${VALID_FILES[@]}"; do
    base="$(basename "$f")"
    module_tag="${base%.genes.regions.tsv}"
    awk -v OFS="\t" -v M="${module_tag}" -F'\t' '
      NF>=2{
        reg=$1; gene=$2
        if (reg ~ /^[^: \t]+:[0-9]+-[0-9]+$/) {
          split(reg,a,/:|-/); s=a[2]+0; e=a[3]+0;
          if (s<1) s=1; if (e<s) next;
          reg=a[1] ":" s "-" e;
          print "coast",M,reg,gene;
          print "north",M,reg,gene;
        } else { bad++ }
      }
      END{ if(bad>0) fprintf(stderr,"[WARN] %d malformed lines in %s (need REGION\\tGENE with REGION like CHR:START-END)\n",bad,FILENAME) }
    ' "$f" >> "${MANIFEST}"
  done
  echo "[OK] Built manifest: ${MANIFEST} with $(wc -l < "${MANIFEST}") rows"
fi

NTASKS=$(wc -l < "${MANIFEST}")
if (( NTASKS == 0 )); then
  echo "[ERR] Empty MANIFEST, nothing to do." >&2
  exit 0
fi

# ========== map array task -> block of MANIFEST rows ==========
TASK=${SLURM_ARRAY_TASK_ID:-0}
if (( TASK < 1 )); then
  echo "[ERR] SLURM_ARRAY_TASK_ID not set" >&2
  exit 1
fi

START=$(( (TASK-1)*ITEMS_PER_TASK + 1 ))
END=$(( TASK*ITEMS_PER_TASK ))
(( END > NTASKS )) && END=$NTASKS

if (( START > NTASKS )); then
  echo "[SKIP] Task ${TASK}: start ${START} > total ${NTASKS} (nothing to do)"
  exit 0
fi

echo "[INFO] Task ${TASK} handling MANIFEST rows ${START}..${END} of ${NTASKS}"

# ==================== ANGSD per MANIFEST row ====================
[[ -r "${REF}" ]] || { echo "[ERR] REF not found: ${REF}" >&2; exit 9; }

# Process each row in this block
sed -n "${START},${END}p" "${MANIFEST}" | while IFS=$'\t' read -r POP MODULE REGION GENE; do
  [[ -n "${POP:-}" && -n "${MODULE:-}" && -n "${REGION:-}" && -n "${GENE:-}" ]] || { echo "[SKIP] malformed row"; continue; }

  # pick bamlist
  if [[ "${POP}" == "coast" ]]; then
    BAMLIST="${COAST_BAMLIST}"
  elif [[ "${POP}" == "north" ]]; then
    BAMLIST="${NORTH_BAMLIST}"
  else
    echo "[WARN] Bad POP in manifest: ${POP}"; continue
  fi
  [[ -s "${BAMLIST}" ]] || { echo "[WARN] Empty BAM list for ${POP}: ${BAMLIST}"; continue; }

  # sanitize filenames
  GENE_TAG="$(printf '%s' "${GENE}" | tr ' /:+\t' '.....')"
  OUTDIR_MOD_POP="${OUTROOT}/${MODULE}.${POP}"
  mkdir -p "${OUTDIR_MOD_POP}"
  OUT="${OUTDIR_MOD_POP}/${GENE_TAG}"

  echo "[INFO] ${POP} | ${MODULE} | ${REGION} | ${GENE}"
  echo "[INFO] OUT=${OUT}"

  # ANGSD
  angsd -b "${BAMLIST}" \
    -ref "${REF}" -anc "${REF}" \
    -r "${REGION}" \
    -doSaf 1 -GL 1 \
    -minMapQ 30 -minQ 20 \
    -P 4 \
    -out "${OUT}"

  # SFS → theta → Tajima D
  if [[ ! -s "${OUT}.saf.idx" ]]; then
    echo "[WARN] no sites for ${GENE}"; continue
  fi

  realSFS "${OUT}.saf.idx" > "${OUT}.sfs" || true
  [[ -s "${OUT}.sfs" ]] || { echo "[WARN] empty SFS for ${GENE}"; continue; }

  realSFS saf2theta "${OUT}.saf.idx" -sfs "${OUT}.sfs" -outname "${OUT}.theta" || true
  [[ -s "${OUT}.theta.thetas.idx" ]] || { echo "[WARN] no thetas for ${GENE}"; continue; }

  thetaStat do_stat "${OUT}.theta.thetas.idx" > "${OUT}.theta.summary.txt"
  thetaStat print   "${OUT}.theta.thetas.idx" > "${OUT}.theta.print.txt"

  echo "[OK] ${MODULE}/${POP} ${GENE} done"
done
