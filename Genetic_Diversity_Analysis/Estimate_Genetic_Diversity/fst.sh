#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END


# Measuring sfs for each population

module load angsd/0.935

angsd -vcf-gl north_pop.vcf -anc ../Ha412HOv2.0-20181130.fasta -dosaf 1 -out north_pop.sfs

angsd -vcf-gl coast_pop.vcf -anc ../Ha412HOv2.0-20181130.fasta -dosaf 1 -out coast_pop.sfs

#calculate all pairwise 2dsfs's


realSFS coast_pop.sfs.saf.idx north_pop.sfs.saf.idx > coastal.north_inland.ml


#---------------------------------------------------------------------------------------------------------------
# prepare the fst for easy analysis etc


realSFS fst index coast_pop.sfs.saf.idx north_pop.sfs.saf.idx -sfs coastal.north_inland.ml -fstout all_pop_pairwise

realSFS fst stats2 pairwise.fst.idx -win 50000 -step 10000 > all_pop_slidingwindow

realSFS fst stats2 pairwise.fst.idx > all_pop_fst

#------------------------------------------------------------------------------------------------------------------
realSFS fst print

realSFS fst print all_pop_pairwise.fst.idx > all_pop_fst_stat



