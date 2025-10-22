# Get NCBI loc IDs and corresponding gff files for each test set

# Extract XP IDs from *_genes.fasta.txt and save to *_ncbi_loc.txt
for file in final_enrichment2/blast_result/*_genes.fasta.txt; do
    # Extract the base name of the file
    base_name="${file##*/}"
    base_name="${base_name%_genes.fasta.txt}"

    # Extract XP IDs and save to a new file
    awk '{print $2}' "$file" | grep -o 'XP_[0-9]*\.[0-9]*' > "final_enrichment2/blast_result/${base_name}_ncbi_loc.txt"
done

# Get gene exon GFF files for selection outlier genes
for file in final_enrichment2/blast_result/*_ncbi_loc.txt; do
    # Extract the base name of the file
    base_name="${file##*/}"
    base_name="${base_name%_ncbi_loc.txt}"

    # Find matching exons and save to a new file
    grep -Ff "$file" all_genes_exons.gff > "final_enrichment2/blast_result/${base_name}_exons.gff"
done
