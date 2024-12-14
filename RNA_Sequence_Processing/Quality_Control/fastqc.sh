#!/bin/bash


## Loading modules required for script commands
module load fastqc/0.11.5

## Running FASTQC
fastqc *.fq*


