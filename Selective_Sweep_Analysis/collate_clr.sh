#!/usr/bin/env bash
#SBATCH -J collate_angsd_pestpg
#SBATCH -c 3                      # one CPU is enough
#SBATCH --mem=2G                  # awk/grep/sed are light; bump if needed
#SBATCH -t 02:00:00               # conservative wall clock
#SBATCH -p cpu                # use a CPU partition on your cluster
#SBATCH -o slurm-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --export=NONE

set -euo pipefail
module purge || true              # clean env (optional)
export LC_ALL=C

# ---- root containing module.pop directories ----
ROOT="/project/pi_brook_moyers_umb_edu/SF2/clr_analysis5/sf2_color_modules_by_gene_out"

# ---- outputs (written in ROOT) ----
SITE_TSV="${ROOT}/sf2_clr_all_sites.tsv"
MAX_TSV="${ROOT}/sf2_clr_max_by_gene.tsv"

# ---- headers ----
echo -e "module\tpop\tgene\tlocation\tLR\talpha\tfile" > "${SITE_TSV}"

# ---- walk module.pop dirs and grab every *_clr.out ----
shopt -s nullglob
for d in "${ROOT}"/*.*; do
  [[ -d "$d" ]] || continue
  modpop="$(basename "$d")"
  module="${modpop%.*}"
  pop="${modpop##*.}"

  # find gene subdirs and their *_clr.out files
  while IFS= read -r -d '' f; do
    gdir="$(dirname "$f")"
    gene="$(basename "$gdir")"
    gene="${gene#gene.}"  # strip "gene." prefix if present

    # append all rows (skip header lines like "location  LR  alpha")
    awk -v MOD="$module" -v POP="$pop" -v GENE="$gene" -v FILE="$f" \
        -v OFS="\t" '
      $0 ~ /^[[:space:]]*$/ { next }                # skip blank
      NR==1 && $1=="location" { next }              # skip header
      NF>=3 {
        loc=$1; lr=$2; alpha=$3;
        print MOD, POP, GENE, loc, lr, alpha, FILE
      }
    ' "$f" >> "${SITE_TSV}"

  done < <(find "$d" -mindepth 2 -maxdepth 2 -type f -name '*_clr.out' -print0)
done

# ---- build per-gene max (max LR + its location) from the site table ----
# Keep one record per module/pop/gene with the greatest LR (tie -> first seen).
awk -F'\t' -v OFS='\t' '
  NR==1 { next }  # skip header of SITE_TSV
  {
    key = $1 FS $2 FS $3               # module\tpop\tgene
    lr = $5 + 0
    if (!(key in best) || lr > best[key]) {
      best[key] = lr
      loc[key]  = $4
      file[key] = $7
    }
  }
  END {
    print "module","pop","gene","max_LR_location","max_LR","file"
    for (k in best) {
      split(k, a, FS)
      print a[1], a[2], a[3], loc[k], best[k], file[k]
    }
  }
' "${SITE_TSV}" \
| sort -t$'\t' -k1,1 -k2,2 -k3,3 > "${MAX_TSV}"

echo "[OK] wrote:"
echo "  ${SITE_TSV}   # all CLR rows (long table)"
echo "  ${MAX_TSV}    # per-gene max LR summary"
