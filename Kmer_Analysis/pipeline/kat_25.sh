#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out kat_k25.log

conda activate kat
kat hist -m 25 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k25_kat.jf
kat gcp -m 25 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k25_kat.gcp
kat comp -m 25 -H 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k25.comp 

kat hist -m 25 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k25_kat.jf
kat gcp -m 25 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k25_kat.gcp
kat comp -m 25 -H 1000000000 -t 10 Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -n -o Friedmanniomyces_simplex_k25_kat.comp
