#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

./NGSadmix -likes genolike.beagle.gz -K 1 -o str_k_1 -P 10

./NGSadmix -likes genolike.beagle.gz -K 2 -o str_k_2 -P 10

./NGSadmix -likes genolike.beagle.gz -K 3 -o str_k_3 -P 10

./NGSadmix -likes genolike.beagle.gz -K 4 -o str_k_4 -P 10

./NGSadmix -likes genolike.beagle.gz -K 5 -o str_k_5 -P 10

./NGSadmix -likes genolike.beagle.gz -K 6 -o str_k_6 -P 10

./NGSadmix -likes genolike.beagle.gz -K 7 -o str_k_7 -P 10
