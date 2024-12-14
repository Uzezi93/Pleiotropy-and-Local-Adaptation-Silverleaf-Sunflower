#!/bin/bash

module load plink/1.07

plink --file argo_plink --recodeAD --out filtered_vcf_plink

