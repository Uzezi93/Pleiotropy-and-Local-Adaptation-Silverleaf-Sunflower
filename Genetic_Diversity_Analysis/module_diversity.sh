#!/bin/bash
#SBATCH -c 5
#SBATCH --mem-per-cpu=8000
#SBATCH -p gpu
#SBATCH -G 4
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out
#SBATCH --mail-type=END

# =============================
# ANGSD Diversity Analysis
# =============================

module load angsd/0.935

GENOME=../../../genome/Ha412HOv2.0-20181130.fasta
OUTDIR=module_diversity/SFS
mkdir -p $OUTDIR

# -----------------------------
# Step 1: SAF for each module
# -----------------------------
for vcf in module_diversity/*.vcf.gz; do
    base=$(basename $vcf .vcf.gz)
    echo "Processing $base"
    angsd -vcf-gl $vcf \
          -anc $GENOME \
          -dosaf 1 \
          -out $OUTDIR/$base.sfs
done

# -----------------------------
# Step 2: Folded SFS
# -----------------------------
for idx in $OUTDIR/*.idx; do
    echo "Folded SFS for $idx"
    realSFS $idx -fold 1 > ${idx}.sfs
done

# -----------------------------
# Step 3: Unfolded SFS
# -----------------------------
for idx in $OUTDIR/*.idx; do
    echo "Unfolded SFS for $idx"
    realSFS $idx > ${idx}_unfold.sfs
done

# -----------------------------
# Step 4: Theta estimation
# -----------------------------
for idx in $OUTDIR/*.idx; do
    echo "Theta estimation for $idx"
    realSFS saf2theta $idx \
        -sfs ${idx}_unfold.sfs \
        -outname $idx
done

# -----------------------------
# Step 5: Print theta stats
# -----------------------------
for thetas in $OUTDIR/*.thetas.idx; do
    echo "Printing theta stats for $thetas"
    thetaStat print $thetas > ${thetas}_theta_stat
done

# -----------------------------
# Step 6: Tajima’s D & others
# -----------------------------
for thetas in $OUTDIR/*.thetas.idx; do
    echo "Running Tajima’s D for $thetas"
    thetaStat do_stat $thetas
done
