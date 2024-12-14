#!/bin/bash

#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END


# perform linkage pruning - i.e. identify prune sites
# plink2 --vcf filtered_snp.recode.vcf.gz \
# --double-id \
# --allow-extra-chr \
# --set-missing-var-ids @:# \
# --indep-pairwise 50 10 0.1 \
# --out argo

# prune and create pca

# Create allele frequency file
./plink2 --bfile data_keep2 --freq --allow-extra-chr --out argo

./plink2 --vcf filtered_snp.recode.vcf.gz \
--double-id \
--allow-extra-chr \
--read-freq argo.afreq \
--set-missing-var-ids @:# \
--make-bed \
--pca \
--out argo
