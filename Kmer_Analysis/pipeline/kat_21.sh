#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out kat_k21.log

conda activate kat
kat hist -m 21 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k21_kat.jf
kat gcp -m 21 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k21.gcp
kat comp -m 21 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k21.comp 

kat hist -m 21 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k21_kat.jf
kat gcp -m 21 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k21.gcp
kat comp -m 21 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -n -o Friedmanniomyces_simplex_k21.comp
