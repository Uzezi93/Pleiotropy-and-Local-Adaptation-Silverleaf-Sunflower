#!/bin/bash

module load vcftools/svn.code.sf.net_r1005

# Convert vcf to plink output

vcftools --vcf filtered_vcf_biallelic.recode.vcf --plink --chrom-map argo.chrom-map.txt --out argo_plink
