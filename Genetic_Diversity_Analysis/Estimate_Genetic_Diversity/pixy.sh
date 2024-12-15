#!/bin/bash

# requires bcftools/bgzip/tabix and vcftools
module load bcftools/2014_06_02
module load vcftools/svn.code.sf.net_r1005
module load tabix/0.2.6
module load pixy/1.2.6

# create a filtered VCF containing only invariant sites
# bgzip output3.vcf
# bgzip filtered_snp.recode.vcf
# tabix output3.vcf.gz

# vcftools --gzvcf output3.vcf.gz \
# --max-maf 0 \
# --recode --stdout | bgzip -c > invariant.vcf.gz

# create a filtered VCF containing only variant sites
# vcftools --gzvcf output3.vcf.gz \
# --mac 1 \
# --recode --stdout | bgzip -c > variant.vcf.gz

# index both vcfs using tabix
# tabix invariant.vcf.gz
# tabix variant.vcf.gz
# tabix filtered_snp.recode.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
filtered_snp.recode.vcf.gz invariant.vcf.gz \
-O z > pixy_filtered.vcf.gz

tabix pixy_filtered.vcf.gz

# run pixy to estimate dxy
source /share/pkg/condas/2018-05-11/bin/activate && conda activate pixy_1.2.6
pixy --stats pi fst dxy \
--vcf pixy_filtered.vcf.gz \
--populations coast_north_pop.txt \
--window_size 1 > dxy_stats
conda deactivate && conda deactivate
