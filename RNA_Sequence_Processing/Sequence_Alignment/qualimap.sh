#!/bin/bash

# load modules
module load qualimap/2.2.1

for i in $(ls *.Aligned.sortedByCoord.out.bam);

do

qualimap bamqc \
-bam $i \
-gff "../helianthus.gtf" \
-outfile "bam_qc" \
-outformat "HTML" 

done
