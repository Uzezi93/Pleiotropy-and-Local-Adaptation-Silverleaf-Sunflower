#!/bin/bash

module load plink/2.00

plink2 --vcf filtered2_snp.recode.vcf --allow-extra-chr --make-bed --out plink_file
