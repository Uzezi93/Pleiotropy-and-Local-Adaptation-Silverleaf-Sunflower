#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

# Measuring sfs for each all samples

# Load angsd module
module load angsd/0.935

# Estimating genetic diversity for all dataset
angsd -vcf-gl all_pop.vcf -anc Ha412HOv2.0-20181130.fasta -dosaf 1 -out SFS2/all_pop.sfs

realSFS SFS2/all_pop.sfs.saf.idx -fold 1 > SFS2/all_pop.sfs

realSFS SFS2/all_pop.sfs.saf.idx  > SFS2/all_pop_unfold.sfs

realSFS saf2theta SFS2/all_pop.sfs.saf.idx -sfs SFS2/all_pop_unfold.sfs -outname SFS2/all_pop

thetaStat print all_pop.thetas.idx  > all_pop_theta_stat

thetaStat do_stat all_pop.thetas.idx
