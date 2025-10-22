
module load bcftools/1.19

# 1) List all sample IDs in the VCF
bcftools query -l filtered_snps.vcf.gz > all.samples

# 2) Make a file with the prefixes of samples I want to be retained (one per line)
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
bcftools view -S keep.samples -Oz -o filtered_snps.subset.vcf.gz filtered_snps.vcf.gz
bcftools index -t filtered_snps.subset.vcf.gz
