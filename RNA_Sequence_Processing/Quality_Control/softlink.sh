#!/bin/bash

for fq in ../raw_fastq/*.fq; do
ln -s $fq
done
