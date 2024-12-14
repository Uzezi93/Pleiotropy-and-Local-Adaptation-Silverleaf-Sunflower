#!/bin/bash

module load star/2.7.6a

for i in $(ls *.trim.fq | sed s/_[12].trim.fq// | sort -u);

do echo;
      
$RUN STAR --genomeDir /project/umb_brook_moyers/argo_network/star_mapping/fastq_soft/star_index \
      --readFilesIn ${i}_1.trim.fq,${i}_2.trim.fq  \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix alignments_STAR/$i.
done


