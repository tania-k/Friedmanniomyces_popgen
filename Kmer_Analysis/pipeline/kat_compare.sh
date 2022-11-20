#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out kat_compare.log

conda activate kat

kat comp -n -t 16 -o compare_Friedmanniomyces 'Friedmanniomyces_endolithicus_?.fastq' 'Friedmanniomyces_simplex_?.fastq' 
