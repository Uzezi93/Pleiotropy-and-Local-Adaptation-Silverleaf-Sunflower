#!/bin/bash

module load trimmomatic/0.39

for infile in arg*
do
outfile=$infile\_trim.fq
java -jar trimmomatic-0.39.jar PE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20
done
