#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load vcftools/0.1.14
module load bcftools/1.15

# Filter north and coast vcf based on positions from LFMM, PCAdapt, eQTL, and eGenes analysis

#vcftools --gzvcf filtered_snp.recode.vcf.gz --positions lfmm3_pos.txt --recode --out lfmm

#vcftools --vcf coast_pop.vcf --positions lfmm3_pos.txt --recode --out coast_lfmm

#vcftools --vcf north_pop.vcf --positions lfmm3_pos.txt --recode --out north_lfmm

#vcftools --gzvcf filtered_snp.recode.vcf.gz --positions pcadapt_pos.txt --recode --out pcadapt

#vcftools --vcf coast_pop.vcf --positions pcadapt_pos.txt --recode --out coast_pcadapt

#vcftools --vcf north_pop.vcf --positions pcadapt_pos.txt --recode --out north_pcadapt

vcftools --gzvcf filtered_snp.recode.vcf.gz --positions eQTL_position2.txt --recode --out eQTL2

vcftools --vcf coast_pop.vcf --positions eQTL_position2.txt --recode --out coast_eQTL2

vcftools --vcf north_pop.vcf --positions eQTL_position2.txt --recode --out north_eQTL2

bcftools view -R eGenes_position2.txt filtered_snp.recode.vcf.gz -o eGenes.vcf

bcftools view -S coast_pop.txt eGenes.vcf > coast_eGenes2.vcf

bcftools view -S north_pop.txt eGenes.vcf > north_eGenes2.vcf

# Get eGenes snp positions
bcftools query -f '%CHROM %POS\n' eGenes.vcf > eGenes_position3.txt

# vcf file for shared loci between PCAdapt and LFMM for CLR analysis
# vcftools --vcf coast_pop.vcf --positions shared_loci_pos.txt --recode --out coast_shared
# vcftools --vcf north_pop.vcf --positions shared_loci_pos.txt --recode --out north_shared

# Create a combined text file with all outlier positions
sort -u lfmm_pos.txt pcadapt_pos.txt eQTL_position2.txt eGenes_position3.txt > combined_positions2.txt

# Turn combined positions text file to tab delimited text file
awk -v OFS='\t' '{$1=$1}1' combined_positions2.txt > combined_positons3.txt

# Remove snp type column
awk '!($3="")' combined_positons3.txt > combined_positons2.txt

# Create control vcf set to confirm outlier test
# Use non-tab delimited file
vcftools --gzvcf filtered_snp.recode.vcf.gz --exclude-positions combined_positions2.txt --recode --recode-INFO-all --out control2

# Divide control vcf into northern and coastal populations
bcftools query -f '%CHROM %POS\n' control2.recode.vcf > control_pos2.txt

# Get control snp positions from control vcf
vcftools --vcf coast_pop.vcf --positions control_pos2.txt --recode --out coast_control2

vcftools --vcf north_pop.vcf --positions control_pos2.txt --recode --out north_control2

# vcftools --remove-indv YOUR_INDIVIDUALS_NAME --vcf your_snps.vcf --recode --out your_filtered_snps.vcf
