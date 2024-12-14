#!/bin/bash

module load vcftools/svn.code.sf.net_r1005

# Calculate allele frequency
# vcftools --gzvcf random_subset.vcf.gz --freq2 --out vcf 

# Calculate mean depth per individual
# vcftools --gzvcf random_subset.vcf.gz --depth --out vcf

# Calculate mean depth per site
# vcftools --gzvcf random_subset.vcf.gz --site-mean-depth --out vcf

# Calculate site quality
# vcftools --gzvcf random_subset.vcf.gz --site-quality --out vcf

# Calculate proportion of missing data per individual
# vcftools --gzvcf random_subset.vcf.gz --missing-indv --out vcf

# Calculate proportion of missing data per site
# vcftools --gzvcf random_subset.vcf.gz --missing-site --out vcf

# Calculate heterozygosity and inbreeding coefficient per individual
# vcftools --gzvcf random_subset.vcf.gz --het --out vcf

# Calculate heterozygosity per site
vcftools --gzvcf filtered2_snp.recode.vcf.recode.vcf --hardy --out vcf_filtered
