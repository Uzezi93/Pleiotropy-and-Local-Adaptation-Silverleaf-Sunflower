#!/bin/bash

module load star/2.7.6a

STAR --runMode genomeGenerate --genomeDir /project/umb_brook_moyers/argo_network/index \
            --genomeFastaFiles /project/umb_brook_moyers/argo_network/index/Ha412HOv2.0-20181130.fasta \
            --sjdbGTFfile /project/umb_brook_moyers/argo_network/star_mapping/fastq_soft/alignments/helianthus.gtf \
            --sjdbOverhang 50 --outFileNamePrefix Helianthus_annuus
