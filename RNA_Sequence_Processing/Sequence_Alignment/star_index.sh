#!/bin/bash
#SBATCH -c 10 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 6:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load star/2.7.11a

STAR --runMode genomeGenerate --genomeDir index \
            --genomeFastaFiles genome/Ha412HOv2.0-20181130.fasta \
            --sjdbGTFfile genome/helianthus.gtf \
            --sjdbOverhang 50 --outFileNamePrefix Helianthus_annuus
