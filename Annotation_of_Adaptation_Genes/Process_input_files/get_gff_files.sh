# Get ncbi gene names (XP IDs) from blast results with Helianthus argophyllus ncbi reference genome
# awk '{print $2}' LFMM_blast_results.txt | grep -o 'XP_[0-9]*\.[0-9]*' > lfmm_ncbi_loc.txt
# awk '{print $2}' pcadapt_blast_results.txt | grep -o 'XP_[0-9]*\.[0-9]*' > pcadapt_ncbi_loc.txt

# Get gene exon gff files for selection outlier genes and WGCNA modules with more than expected LFMM outlier genes
# grep -Ff lfmm_ncbi_loc.txt all_genes_exons.gff > lfmm_exons.gff
# grep -Ff pcadapt_ncbi_loc.txt all_genes_exons.gff > pcadapt_exons.gff


# Extract XP IDs from *_genes.fasta.txt and save to *_ncbi_loc.txt
for file in final_enrichment/blast_result/*_genes.fasta.txt; do
    # Extract the base name of the file
    base_name="${file##*/}"
    base_name="${base_name%_genes.fasta.txt}"

    # Extract XP IDs and save to a new file
    awk '{print $2}' "$file" | grep -o 'XP_[0-9]*\.[0-9]*' > "final_enrichment/blast_result/${base_name}_ncbi_loc.txt"
done

# Get gene exon GFF files for selection outlier genes
for file in final_enrichment/blast_result/*_ncbi_loc.txt; do
    # Extract the base name of the file
    base_name="${file##*/}"
    base_name="${base_name%_ncbi_loc.txt}"

    # Find matching exons and save to a new file
    grep -Ff "$file" all_genes_exons.gff > "final_enrichment/blast_result/${base_name}_exons.gff"
done
