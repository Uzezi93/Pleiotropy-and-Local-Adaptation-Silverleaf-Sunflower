#!/bin/bash

module load vcftools/svn.code.sf.net_r1005

# Get raw snp files

# vcftools --gzvcf output.vcf.gz --remove-indels --recode --recode-INFO-all --out raw_snps.vcf.gz

# Add vcf filters

# vcftools --gzvcf raw_snps.vcf.gz.recode.vcf.gz --minQ 30 --minGQ 20 --maf 0.05 --max-non-ref-af 0.999 --hwe 0.8 --max-missing 0.9 --min-meanDP 5 --max-meanDP 40 --recode --recode-INFO-all --out filtered_snp


# Filter north and coast vcf based on positions from LFMM, PCAdapt, eQTL, and eGenes analysis

# vcftools --gzvcf annoted_snps.vcf.gz --positions lfmm_pos.txt --recode --out lfmm

# vcftools --vcf coast_pop.vcf --positions lfmm_pos.txt --recode --out coast_lfmm

# vcftools --gzvcf annoted_snps.vcf.gz --positions pcadapt_pos.txt --recode --out pcadapt

# vcftools --vcf coast_pop.vcf --positions pcadapt_pos.txt --recode --out coast_pcadapt

# vcftools --vcf annoted_snps.vcf --positions eQTL_position.txt --recode --out eQTL

# vcftools --vcf coast_pop.vcf --positions eQTL_position.txt --recode --out coast_eQTL

# vcftools --vcf north_pop.vcf --positions eQTL_position.txt --recode --out north_eQTL

# vcftools --gzvcf annoted_snps.vcf.gz --bed eGenes.bed --out eGenes --recode --keep-INFO-all

# bcftools view -S coast_pop.txt eGenes.recode.vcf > coast_eGenes.vcf

# bcftools view -S north_pop.txt eGenes.recode.vcf > north_eGenes.vcf



# Create control vcf set to confirm outlier test
# vcftools --vcf annoted_snps.vcf --exclude-positions combined_positions.txt --recode --recode-INFO-all --out control

# Create 012 coded vcf file for eQTL analysis
# vcftools --vcf filtered_snp.recode.vcf --012 --recode-INFO-all --out filtered_vcf


# Make CLR input files from LFMM, PCAdapt, eQTL, and eGenes loci from north and coastal populations

cd ./outputs/

 for i in *.vcf

 do

 vcftools --vcf $i --counts2 --out ${i}_SF_tmp

 tail -n+2 ${i}_SF_tmp.frq.count | awk -v OFS="\t" '{print $2,$6,$4,"1"}'> ${i}_SF.in

 echo -e 'position\tx\tn\tfolded' | cat - ${i}_SF.in > temp && mv temp ${i}_SF.in

 done




