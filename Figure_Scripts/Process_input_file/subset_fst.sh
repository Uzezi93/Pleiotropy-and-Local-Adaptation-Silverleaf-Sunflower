# Change position files for outliers and control genes to tab file before subsetting fst

awk -v OFS='\t' '{$1=$1}1' control_pos2.txt > control_pos3.txt

awk -v OFS='\t' '{$1=$1}1' eGenes_position3.txt > eGenes_position4.txt

awk -v OFS='\t' '{$1=$1}1' eQTL_position2.txt > eQTL_position3.txt

# Remove snp column from lfmm postions txt file

#awk '!($3="")' lfmm3_pos.txt > lfmm_pos2.txt

#awk -v OFS='\t' '{$1=$1}1' lfmm_pos2.txt > lfmm_pos.txt

#awk '!($3="")' pcadapt_pos.txt > pcadapt_pos2.txt

#awk -v OFS='\t' '{$1=$1}1' pcadapt_pos2.txt > pcadapt_pos.txt


# Subset fst file into control, eGenes, eQTL, LFMM and PCAdapt outliers

grep -f control_pos3.txt fst_stat > control_fst.txt

grep -f eGenes_position4.txt fst_stat > eGenes_fst.txt

grep -f eQTL_position3.txt fst_stat > eQTL_fst.txt

#grep -f lfmm_pos.txt fst_stat > lfmm_fst.txt

#grep -f pcadapt_pos.txt fst_stat > pcadapt_fst.txt
