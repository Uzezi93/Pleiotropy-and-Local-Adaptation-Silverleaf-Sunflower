#!/bin/bash

module load freebayes/1.3.1

source /share/pkg/condas/2018-05-11/bin/activate && conda activate freebayes_1.3.1

ref="/project/umb_brook_moyers/argo_network/index/Ha412HOv2.0-20181130.fasta"

freebayes-parallel \
   <(/share/pkg/condas/2018-05-11/envs/freebayes_1.3.1/bin/fasta_generate_regions.py ${ref}.fai 100000) 20 \
   --fasta-reference ${ref} \
   --report-monomorphic \
   --bam-list bam.fofn  > output_2.vcf

conda deactivate && conda deactivate
