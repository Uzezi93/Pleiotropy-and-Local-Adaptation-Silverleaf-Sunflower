#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

# We ran this script to compute Fay'H per gene
# Measuring sfs for each population

module load angsd/0.935

angsd -vcf-gl north_lfmm.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_lfmm.sfs

angsd -vcf-gl coast_lfmm.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_lfmm.sfs

angsd -vcf-gl north_pcadapt.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_pcadapt.sfs

angsd -vcf-gl coast_pcadapt.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_pcadapt.sfs

angsd -vcf-gl coast_eQTL.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_eQTL.sfs

angsd -vcf-gl north_eQTL.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_eQTL.sfs

angsd -vcf-gl coast_eGenes.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_eGenes.sfs

angsd -vcf-gl north_eGenes.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_eGenes.sfs

angsd -vcf-gl coast_control.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_control.sfs

angsd -vcf-gl north_control.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_control.sfs

angsd -vcf-gl coast_shared.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/coast_shared.sfs

angsd -vcf-gl north_shared.recode.vcf -anc ../../../genome/Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS/north_shared.sfs
#---------------------------------------------------------------------------------------------------------

# Measuring folded saf for each population

 for i in SFS/*.idx
 do
  realSFS $i -fold 1 > $i.sfs
  done

#----------------------------------------------------------------------------------------------------------
# Measuring unfolded saf for each population

 for i in SFS/*.idx
 do
 realSFS $i > ${i}_unfold.sfs
 done

#----------------------------------------------------------------------------------------------------------
# Calculate theta for each site

 for i in SFS/*.idx
 do
 realSFS saf2theta $i -sfs ${i}_unfold.sfs -outname $i
 done

#----------------------------------------------------------------------------------------------------------
# Print theta statistics

 for i in SFS/*.thetas.idx
 do
 thetaStat print $i > ${i}_theta_stat
 done

#----------------------------------------------------------------------------------------------------------
# Estimate Tajimas D and other statistics

 for i in SFS/*.thetas.idx
 do
 thetaStat do_stat $i
 done

#----------------------------------------------------------------------
