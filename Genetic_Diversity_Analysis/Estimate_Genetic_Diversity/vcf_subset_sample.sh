#!/bin/bash
#SBATCH -J vcf_merge
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 12:00:00
#SBATCH -p gpu                   # if you have a CPU partition, prefer that
#SBATCH -o slurm-%j.out
#SBATCH -G 4                   # GPUs not needed

# Script to filter vcf file with monomorphic sites to contain my specific samples

module load bcftools/1.19

# 1) List all sample IDs in the VCF
bcftools query -l freebayes.merged.norm.vcf.gz > all.samples

# 2) Make a file with the prefixes you want (one per line)
cat > prefixes.txt <<'EOF'
btm9
btm10
btm13
arg6B
btm17
btm19
btm20
btm21
btm22
btm25
btm26
btm27
arg4B
arg14B
btm30
btm31
arg2B
btm32
btm34
EOF

# 3) Keep any sample whose prefix (before first '-') is in that list
awk 'NR==FNR { keep[$1]=1; next }
     { base=$1; sub(/-.*/,"",base); if (keep[base]) print $1 }' \
    prefixes.txt all.samples > keep.samples

# (Sanity check: print the kept sample IDs)
cat keep.samples

# 4) Subset the VCF and index
bcftools view -S keep.samples -Oz -o all_filtered_snps.subset.vcf.gz freebayes.merged.norm.vcf.gz
bcftools index -t all_filtered_snps.subset.vcf.gz
