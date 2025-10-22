#!/usr/bin/env bash
#SBATCH -J pixy_genes
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 48:00:00
#SBATCH -p gpu
#SBATCH -G 12
#SBATCH -o slurm-%A_%a.out
#SBATCH --array=1-1000%24

set -euo pipefail

# -------- config --------
VCF="pixy_diploidized.vcf.gz"      # diploidized, bgzip+tabix indexed
POPS="keep.samples.txt"            # 2-column, TAB-separated: sample <TAB> population
BED="gene_coord.bed.gz"            # 0-based BED; col4 = gene ID
OUTDIR="./pixy_out_genes"
PREFIX_BASE="pixy_genes"
CONTIGS_FILE="contigs.txt"
CHUNK_SIZE=30                      # each array task processes 30 contigs

mkdir -p "${OUTDIR}"

# ---- sanity: VCF must be bgzip + indexed ----
if ! bgzip -t "${VCF}" >/dev/null 2>&1; then
  echo "[ERR] ${VCF} is not valid BGZF (or truncated)." >&2
  exit 2
fi
[[ -f "${VCF}.tbi" || -f "${VCF%.gz}.csi" ]] || bcftools index -ft "${VCF}"

# ---- sanity: BED must be tabix-indexed (BED/TSV index) ----
if [[ ! -f "${BED}.tbi" ]]; then
  echo "[INFO] Indexing BED..."
  tabix -f -p bed "${BED}"
fi

# Build contig list from VCF if missing
if [[ ! -s "${CONTIGS_FILE}" ]]; then
  bcftools view -h "${VCF}" | awk -F'[=,]' '/^##contig=/{print $3}' > "${CONTIGS_FILE}"
fi
N_CONTIGS=$(wc -l < "${CONTIGS_FILE}")
(( N_CONTIGS > 0 )) || { echo "[ERR] No contigs found in ${VCF}"; exit 3; }

# ---- compute this array task's contig block ----
TASK=${SLURM_ARRAY_TASK_ID:-1}
START=$(( (TASK-1)*CHUNK_SIZE + 1 ))
END=$(( TASK*CHUNK_SIZE ))
(( START <= N_CONTIGS )) || { echo "[SKIP] Task ${TASK}: start ${START} > ${N_CONTIGS}"; exit 0; }
(( END > N_CONTIGS )) && END=${N_CONTIGS}

echo "[pixy-array] Task ${TASK}: contigs ${START}-${END} of ${N_CONTIGS}"

# ---- loop contigs ----
for ((i=START; i<=END; i++)); do
  CONTIG=$(sed -n "${i}p" "${CONTIGS_FILE}")
  echo "[pixy] ${TASK}: ${CONTIG}"

  # Extract only this contig's gene intervals to a temp BED (skip if none)
  TMPBED=$(mktemp)
  # tabix returns headerless BED lines; we don't need a header
  if ! tabix "${BED}" "${CONTIG}" > "${TMPBED}" 2>/dev/null; then
    echo "[WARN] No BED intervals for ${CONTIG}; skipping."
    rm -f "${TMPBED}"
    continue
  fi
  if [[ ! -s "${TMPBED}" ]]; then
    echo "[WARN] Empty BED for ${CONTIG}; skipping."
    rm -f "${TMPBED}"
    continue
  fi

  # Run pixy per gene on this contig (no --window_size when using --bed_file)
  pixy --stats pi dxy fst watterson_theta tajima_d \
       --vcf "${VCF}" \
       --populations "${POPS}" \
       --chromosomes "${CONTIG}" \
       --bed_file "${TMPBED}" \
       --output_prefix "${PREFIX_BASE}_${CONTIG}" \
       --output_folder "${OUTDIR}"

  rm -f "${TMPBED}"
done

echo "[pixy-array] Task ${TASK} done."
