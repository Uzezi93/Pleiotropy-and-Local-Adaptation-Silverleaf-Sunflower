#!/bin/bash

# Generate read counts for Deseq and WGCNA
module load subread/1.6.2

mkdir Counts
featureCounts -T 20 -p -B -O --fraction -a helianthus.gtf -f -o counts.txt *.rem_dups_w_mate.bam

