#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out simplex_gcp.log

kat gcp -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex.gcp
