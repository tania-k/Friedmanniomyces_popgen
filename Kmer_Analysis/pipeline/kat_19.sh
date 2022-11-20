#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out kat_k19.log

conda activate kat
kat hist -m 19 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k19_kat.jf
kat gcp -m 19 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k19.gcp
kat comp -m 19 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k19.comp 

kat hist -m 19 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k19_kat.jf
kat gcp -m 19 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k19.gcp
kat comp -m 19 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -n -o Friedmanniomyces_simplex_k19.comp
