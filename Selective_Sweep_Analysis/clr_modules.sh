#!/usr/bin/env bash
#SBATCH -J sf2_color_modules_by_gene
#SBATCH -c 5
#SBATCH --mem-per-cpu=8000
#SBATCH -t 10:00:00
#SBATCH -p gpu
#SBATCH -G 4
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=END
#SBATCH --array=1-1000%1000

set -euo pipefail

# ==================== batching config ====================
# How many MANIFEST rows each array task should process?
ITEMS_PER_TASK=60

# -------------------- modules --------------------
module load bcftools/1.19
#module load vcftools/0.1.13

# -------------------- PATHS --------------------
# Color modules root
ROOT="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/enriched_modules_gene_lists"

# per-pop VCFs named <module>.<pop>.vcf.gz (coast|north)
VCF_DIR="${ROOT}/by_pop_vcfs"

# per-module BEDs named <module>.coords.bed
BED_DIR="${ROOT}"

# Output root and manifest
OUTROOT="${PWD}/sf2_color_modules_by_gene_out"
MANIFEST="${OUTROOT}/task_manifest.tsv"

# SweepFinder2 binary
SF2="/project/pi_brook_moyers_umb_edu/SF2/SweepFinder2"

mkdir -p "${OUTROOT}"

# -------------------- MODULE & BED MAPPING --------------------
# match the VCF prefix "<module>.<pop>.vcf.gz"
# and the BED filenames listed.
declare -A BED_FOR=(
  [blue_genes]="${BED_DIR}/blue_genes.coords.bed"
  [green_genes]="${BED_DIR}/green_genes.coords.bed"
  [pink_genes]="${BED_DIR}/pink_genes.coords.bed"
  [red_genes]="${BED_DIR}/red_genes.coords.bed"
  [yellow_genes]="${BED_DIR}/yellow_genes.coords.bed"
  [control_pool_all_genes]="${BED_DIR}/control_pool_all_genes.coords.bed"
)

# modules and populations to process
MODULES=(blue_genes green_genes pink_genes red_genes yellow_genes control_pool_all_genes)
POPS=(coast north)

# -------------------- BUILD MANIFEST --------------------
# MANIFEST columns: module  pop  chrom  start  end  gene  vcf
if [[ ! -s "${MANIFEST}" ]]; then
  echo -e "module\tpop\tchrom\tstart\tend\tgene\tvcf" > "${MANIFEST}.tmp"

  for mod in "${MODULES[@]}"; do
    BED="${BED_FOR[$mod]:-}"
    if [[ -z "${BED}" || ! -s "${BED}" ]]; then
      echo "[WARN] Missing BED for ${mod}: ${BED}" >&2
      continue
    fi
    sed -i 's/\r$//' "${BED}"  # normalize line endings

    for pop in "${POPS[@]}"; do
      VCF="${VCF_DIR}/${mod}.${pop}.vcf.gz"
      if [[ ! -s "${VCF}" ]]; then
        echo "[WARN] Missing VCF for ${mod}/${pop}: ${VCF}" >&2
        continue
      fi
      # BED with at least 4 columns: chrom start end gene (extra columns okay)
      awk -v M="${mod}" -v P="${pop}" -v V="${VCF}" -v OFS="\t" '
        NF>=4 {print M, P, $1, $2, $3, $4, V}
      ' "${BED}" >> "${MANIFEST}.tmp"
    done
  done

  mv "${MANIFEST}.tmp" "${MANIFEST}"
  echo "[OK] Built manifest: ${MANIFEST} with $(($(wc -l < "${MANIFEST}")-1)) tasks"
fi

# -------------------- SLICE this array task --------------------
TOTAL_ROWS=$(( $(wc -l < "${MANIFEST}") - 1 ))   # minus header
IDX=${SLURM_ARRAY_TASK_ID:-0}
if (( IDX < 1 )); then
  echo "[ERR] Bad SLURM_ARRAY_TASK_ID=${IDX}"; exit 2
fi

SLICE_START=$(( (IDX - 1) * ITEMS_PER_TASK + 1 ))  # 1-based over data rows
SLICE_END=$(( SLICE_START + ITEMS_PER_TASK - 1 ))
(( SLICE_END > TOTAL_ROWS )) && SLICE_END=${TOTAL_ROWS}

if (( SLICE_START > TOTAL_ROWS )); then
  echo "[SKIP] Array index ${IDX} has no work (START=${SLICE_START} > TOTAL_ROWS=${TOTAL_ROWS})"
  exit 0
fi

echo "[INFO] Array ${IDX}: processing manifest rows ${SLICE_START}..${SLICE_END} of ${TOTAL_ROWS}"

# -------------------- per-gene runner --------------------
run_one_gene () {
  local MODULE="$1" POP="$2" CHR="$3" GSTART="$4" GEND="$5" GENE="$6" VCF="$7"

  # sanity
  [[ -s "${VCF}" ]] || { echo "[WARN] missing VCF ${VCF}; skipping ${GENE}"; return 0; }

  # sanitize gene id for filenames
  local GENE_TAG; GENE_TAG="$(printf '%s' "${GENE}" | tr ' /:+\t' '.....')"
  local WORK="${OUTROOT}/${MODULE}.${POP}/${GENE_TAG}"
  mkdir -p "${WORK}"

  echo "[INFO] ${MODULE}/${POP}  ${CHR}:${GSTART}-${GEND}  ${GENE}"

  # subset vcf
  local SUBVCF="${WORK}/${GENE_TAG}.vcf.gz"
  bcftools view -r "${CHR}:${GSTART}-${GEND}" -Oz -o "${SUBVCF}" "${VCF}" 2> "${WORK}/bcftools.view.log" || true
  [[ -s "${SUBVCF}" ]] && bcftools index -f "${SUBVCF}" || true

  # variants check
  local NVAR; NVAR=$( [[ -s "${SUBVCF}" ]] && bcftools index -n "${SUBVCF}" || echo 0 )
  if (( NVAR == 0 )); then
    echo "[WARN] no variants in ${GENE}; skipping"
    return 0
  fi

  # counts2 -> .in
  local TMPPREFIX="${WORK}/${GENE_TAG}_SF_tmp"
  vcftools --gzvcf "${SUBVCF}" --counts2 --out "${TMPPREFIX}" > "${WORK}/vcftools.counts2.log" 2>&1 || true
  local FRQ="${TMPPREFIX}.frq.count"
  if [[ ! -s "${FRQ}" ]]; then
    echo "[WARN] empty frq.count for ${GENE}; skipping"
    return 0
  fi

  local INFILE="${WORK}/${GENE_TAG}.in"
  {
    echo -e "position\tx\tn\tfolded"
    # POS=col2, n=col4, x=col6, folded=1
    tail -n +2 "${FRQ}" | awk -v OFS='\t' '{print $2,$6,$4,"1"}'
  } > "${INFILE}"

  # grid
  local GRID="${WORK}/${GENE_TAG}_grid.txt"
  bcftools query -f '%POS\n' "${SUBVCF}" > "${GRID}"

  # SweepFinder2
  local SPEC="${WORK}/${GENE_TAG}_specfile.out"
  local CLR="${WORK}/${GENE_TAG}_clr.out"
  "${SF2}" -f "${INFILE}" "${SPEC}" > "${WORK}/sf2_spec.log" 2>&1 || true
  "${SF2}" -lu "${GRID}" "${INFILE}" "${SPEC}" "${CLR}" > "${WORK}/sf2_run.log" 2>&1 || true

  echo "[OK] ${MODULE}/${POP} ${GENE} -> ${CLR}"
}

# -------------------- iterate over this task's slice --------------------
awk -v s="${SLICE_START}" -v e="${SLICE_END}" 'NR==1{next} {i=NR-1; if(i>=s && i<=e) print $0}' "${MANIFEST}" \
| while IFS=$'\t' read -r MODULE POP CHR START END GENE VCF; do
    [[ -z "${MODULE}" || -z "${VCF}" ]] && { echo "[WARN] bad row; skipping"; continue; }
    run_one_gene "${MODULE}" "${POP}" "${CHR}" "${START}" "${END}" "${GENE}" "${VCF}"
  done
