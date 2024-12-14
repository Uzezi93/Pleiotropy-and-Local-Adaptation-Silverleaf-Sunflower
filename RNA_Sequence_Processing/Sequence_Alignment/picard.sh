#!/bin/bash

module load picard/2.23.3

#for i in $(ls *.Aligned.sortedByCoord.out.bam);

#do

#java -jar picard.jar AddOrReplaceReadGroups \
 #     I=${i} \
  #    O=${i}-sorted-md-rg.bam \
   #   RGID=${i} \
    #  RGLB=${i} \
     # RGPL=illumina \
      #RGPU=unit1 \
      #RGSM=${i}

#done

#for i in $(ls *.Aligned.sortedByCoord.out.bam-sorted-md-rg.bam.rem_dups_w_mate.bam.g.vcf);

 #do

java -jar picard.jar VcfToIntervalList \
      I=btm9-4.Aligned.sortedByCoord.out.bam-sorted-md-rg.bam.rem_dups_w_mate.bam.g.vcf 
      O=btm9.interval_list
#done
