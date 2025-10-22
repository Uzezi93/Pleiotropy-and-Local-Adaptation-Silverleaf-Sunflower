#!/bin/bash
#SBATCH -J pixy_vcf_filter
#SBATCH -c 2
#SBATCH --mem=80G
#SBATCH -t 24:00:00
#SBATCH -p cpu
#SBATCH -o slurm-%j.out

module load bcftools/1.19

# 1) Diploidize GTs everywhere
( bcftools view -h pixy_filtered.vcf.gz
  bcftools view -H pixy_filtered.vcf.gz \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        # rewrite per-sample fields
        for(i=10;i<=NF;i++){
          n = split($i, a, ":"); gt=a[1]
          if     (gt == ".")          a[1] = "./."
          else if (gt ~ /^[0-9]+$/)   a[1] = gt "/" gt
          # else: already diploid (0/0, 0/1, 1|1, ./., .|.)
          $i = a[1]; for(j=2;j<=n;j++) $i = $i ":" a[j]
        }
        print
      }'
) | bgzip > pixy_diploidized.vcf.gz
bcftools index -t -f  pixy_diploidized.vcf.gz

# 2) keep sites with >=2 called genotypes
bcftools view -i 'N_MISSING<=NSAMPLES-2' -Oz -o pixy_diploidized.min2.vcf.gz pixy_diploidized.vcf.gz
bcftools index -t -f pixy_diploidized.min2.vcf.gz
