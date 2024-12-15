#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G 4 # Number of GPUs
#SBATCH -t 10:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load bcftools/1.15
module load vcftools/0.1.13

# Subset vcf into scaffolds
bcftools query -f '%CHROM\n' north_lfmm.recode.vcf.gz | uniq > saffolds.txt

cat saffolds.txt | while read CONTIG; do
    vcftools --gzvcf north_lfmm.recode.vcf.gz  --chr $CONTIG --recode --recode-INFO-all --out ./scaff/$CONTIG.out
done

cd scaff

# make clr input files

for i in *.vcf; do
  # Run vcftools to generate counts
  vcftools --vcf "${i}" --counts2 --out "${i%.vcf}_tmp"

  # Process the output file, tailing from the second line, and format with awk
  tail -n +2 "${i%.vcf}_tmp.frq.count" | awk -v OFS='\t' '{print $2, $6, $4, "1"}' > "${i%.vcf}.input"

  # Add header and combine with the processed input
  echo -e 'position\tx\tn\tfolded' | cat - "${i%.vcf}.input" > temp && mv temp "${i%.vcf}.input"
done

# Get grid files for CLR analysis

for i in *.vcf; do
  bcftools query -f '%POS\n' "${i}" > "${i%.vcf}_grid.txt"
done


# Repair grid files

#for i in *.txt; do

# sed 's/^.\{4\}//g' "${i}" > "${i%.txt}.txt"; done


# Calculate CLR

# Create specfile

for i in *.input; do
/project/pi_brook_moyers_umb_edu/SF2/SweepFinder2 -f "${i}" "${i%.input}_specfile.out"; done


# Loop through each input file with the specified suffixes
for grid_file in *_grid.txt; do
    # Extract the filename without the directory and suffix
    filename=$(basename "$grid_file" _grid.txt)

    # Form the names of the other two files
    specfile="${filename}_specfile.out"
    inputfile="${filename}.input"

    # Run SweepFinder2 with the three files
    /project/pi_brook_moyers_umb_edu/SF2/SweepFinder2 -lu "$grid_file" "$inputfile" "$specfile" "${filename}_clr.out"

done

# combine all clr outputs to one file
# head -n 1 Ha412HOChr17.out.recode_clr.out  > coast_control_clr.out

# tail -n +2 -q *clr.out >> coast_control_clr.out

cat *_clr.out > temp_clr_combined.txt

awk '/location[[:space:]]+LR[[:space:]]+alpha/ {if (seen++) next} {print}' temp_clr_combined.txt > coast_control2_clr.out
