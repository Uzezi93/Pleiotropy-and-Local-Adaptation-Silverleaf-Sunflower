#!/bin/bash
#SBATCH -c 5 # Number of Cores per Task
#SBATCH --mem-per-cpu=8000  # Requested Memory
#SBATCH -p gpu  # Partition
#SBATCH -G  4 # Number of GPUs
#SBATCH -t 01:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = Job ID
#SBATCH --mail-type=END

module load ncbi-rmblastn/2.11.0
module load seqtk/1.3

# Get fasta file for all genes in our expression data
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta all_genes_expressed.txt >  all_genes.fasta

# Get fasta file for LFMM genes
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta LFMM_genes.txt >  LFMM_genes.fasta

# Get fasta file for pcadapt genes
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta pcadapt_genes.txt >  pcadapt_genes.fasta

# Get fasta file for genes in ivory module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta ivory_genes.txt >  ivory_genes.fasta

# Get fasta file for genes in steelblue module
# seqtk subseq Ha412HOv2.0-20181130.transcripts.fasta steelblue_genes.txt >  steelblue_genes.fasta

# Create a BLAST database from the sunflower genome FASTA file
# makeblastdb -in NCBI/GCF_002127325.2_HanXRQr2.0-SUNRISE_cds_from_genomic.fa -dbtype nucl -out sunflower_db

# Blast the all_genes.fasta query against the sunflower genome database
# blastn -query all_genes.fasta -db sunflower_db -outfmt "6 qseqid sacc" -out all_genes_blast_results.txt

# Blast the LFMM_genes.fasta query against the sunflower genome database
# blastn -query LFMM_genes.fasta -db sunflower_db -outfmt "6 qseqid sacc" -out LFMM_blast_results.txt
# blastn -query pcadapt_genes.fasta -db sunflower_db -outfmt "6 qseqid sacc" -out pcadapt_blast_results.txt

# Blast module genes fasta file against sunflower genome database
# blastn -query ivory_genes.fasta -db sunflower_db -outfmt "6 qseqid sacc" -out ivory_blast_results.txt
# blastn -query steelblue_genes.fasta -db sunflower_db -outfmt "6 qseqid sacc" -out steelblue_blast_results.txt

# Blast the module fasta files  query against the sunflower genome database
 for i in modules/fasta_files/*.fasta; do
   blastn -query "$i" -db sunflower_db -outfmt "6 qseqid sacc" -out modules/blast_results/ "$i.txt"
done