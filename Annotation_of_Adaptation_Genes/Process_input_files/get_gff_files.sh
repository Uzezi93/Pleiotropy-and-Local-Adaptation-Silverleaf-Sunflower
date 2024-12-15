# Get ncbi gene names (XP IDs) from blast results with Helianthus argophyllus ncbi reference genome
awk '{print $2}' LFMM_blast_results.txt | grep -o 'XP_[0-9]*\.[0-9]*' > lfmm_ncbi_loc.txt
awk '{print $2}' pcadapt_blast_results.txt | grep -o 'XP_[0-9]*\.[0-9]*' > pcadapt_ncbi_loc.txt

# Extract NCBI gene names (XP IDs) for each module file
for file in modules/ *_blast_results.txt; do
    # Extract the base name of the file
    base_name=$(basename "$file" "_blast_results.txt")

    # Extract XP IDs and save to a new file
    awk '{print $2}' "$file" | grep -o 'XP_[0-9]*\.[0-9]*' > "modules/${base_name}_ncbi_loc.txt"
done

# Get gene exon gff files for selection outlier genes and WGCNA modules with more than expected outlier genes
grep -Ff lfmm_ncbi_loc.txt all_genes_exons.gff > lfmm_exons.gff
grep -Ff pcadapt_ncbi_loc.txt all_genes_exons.gff > pcadapt_exons.gff

for file in modules/ *_ncbi_loc.txt; do
    # Extract the base name of the file
    base_name=$(basename "$file" "_ncbi_loc.txt")

    # Extract XP IDs and save to a new file
    awk '{print $2}' "$file" | grep -o 'XP_[0-9]*\.[0-9]*' > "modules/${base_name}_exons.gff"
done
