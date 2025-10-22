#!/usr/bin/env bash
#SBATCH -J make_clr_inputs_by_pop
#SBATCH -c 5
#SBATCH --mem-per-cpu=8000
#SBATCH -t 01:30:00
#SBATCH -p gpu
#SBATCH -G 4
#SBATCH -o slurm-%A_%a.out
#SBATCH --mail-type=END
#SBATCH --array=1-6     # override at submit: sbatch --array=1-N%X this_script.sh

set -euo pipefail

module load vcftools/0.1.16
module load htslib/1.19.1 || true   # bgzip/tabix
module load bcftools/1.19   || true

# -------- paths --------
INDIR="/project/pi_brook_moyers_umb_edu/Uzezi_argo/raw_fastq/variants/chunks/analysis_vcf_files/vcf_by_gene_sets"
OUTVCF="${INDIR}/by_pop_vcfs"         # subset VCFs per module×pop
OUTSF="${INDIR}/clr_inputs"           # CLR/SweepFinder .in files per module×pop
COAST="${INDIR}/../coast.samples"     # adjust if stored elsewhere
NORTH="${INDIR}/../north.samples"

mkdir -p "${OUTVCF}" "${OUTSF}"

# -------- sanity --------
[[ -s "${COAST}" ]] || { echo "[ERR] coast.samples not found: ${COAST}" >&2; exit 2; }
[[ -s "${NORTH}" ]] || { echo "[ERR] north.samples not found: ${NORTH}" >&2; exit 3; }

shopt -s nullglob
# Build a deterministic list of VCFs
mapfile -t vcfs < <(printf "%s\n" "${INDIR}/"*.vcf.gz | sort -V)
NV=${#vcfs[@]}
if (( NV == 0 )); then
  echo "[ERR] No VCFs matched in ${INDIR}" >&2
  exit 4
fi

IDX=${SLURM_ARRAY_TASK_ID:-0}
if (( IDX < 1 || IDX > NV )); then
  echo "[SKIP] SLURM_ARRAY_TASK_ID=${IDX} outside 1..${NV}"
  exit 0
fi

vcf="${vcfs[$((IDX-1))]}"
base=$(basename "${vcf}" .vcf.gz)   # e.g., control_genes
echo "=========================================================="
echo "[INFO] Task ${IDX}/${NV} processing: ${base}"

# -------- helper: make subset vcf --------
make_subset_vcf () {
  local vcf_in="$1"
  local keep_samples="$2"
  local out_stem="$3"   # without .vcf.gz

  echo "[INFO] subset VCF -> ${OUTVCF}/${out_stem}.vcf.gz"
  vcftools \
    --gzvcf "${vcf_in}" \
    --keep "${keep_samples}" \
    --remove-indels \
    --min-alleles 2 --max-alleles 2 \
    --recode --recode-INFO-all --stdout \
  | bgzip -c > "${OUTVCF}/${out_stem}.vcf.gz"

  tabix -f -p vcf "${OUTVCF}/${out_stem}.vcf.gz" || true
}

# -------- helper: counts2 -> CLR .in --------
make_clr_in () {
  local vcf_path="$1"
  local out_stem="$2"   # without .in

  echo "[INFO] counts2 -> ${OUTSF}/${out_stem}.frq.count"
  vcftools \
    --gzvcf "${vcf_path}" \
    --counts2 \
    --out "${OUTSF}/${out_stem}_SF_tmp"

  local frq="${OUTSF}/${out_stem}_SF_tmp.frq.count"
  if [[ ! -s "${frq}" ]]; then
    echo "[WARN] Empty counts for ${out_stem}; skipping .in"
    return 0
  fi

  # CLR/SweepFinder (folded): position  x  n  folded
  echo -e "position\tx\tn\tfolded" > "${OUTSF}/${out_stem}.in"
  tail -n +2 "${frq}" \
    | awk -v OFS="\t" '{print $2,$6,$4,"1"}' \
    >> "${OUTSF}/${out_stem}.in"

  local nlines
  nlines=$(($(wc -l < "${OUTSF}/${out_stem}.in") - 1))
  echo "[OK] ${out_stem}.in rows (excluding header): ${nlines}"
}

# -------- do the work for this one VCF --------
make_subset_vcf "${vcf}" "${COAST}" "${base}.coast"
make_clr_in     "${OUTVCF}/${base}.coast.vcf.gz" "${base}.coast"

make_subset_vcf "${vcf}" "${NORTH}" "${base}.north"
make_clr_in     "${OUTVCF}/${base}.north.vcf.gz" "${base}.north"

echo "[DONE] ${base}"
