#!/bin/bash

# module load plink/2.00
module load plink/1.07

# convert vcf file to binary file for eQTL analysis
plink --file argo_plink --recode12 --out eQTL_file2
