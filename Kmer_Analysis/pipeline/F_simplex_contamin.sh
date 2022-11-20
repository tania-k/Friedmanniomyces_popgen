#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out contam_simp_k21.log

kat sect -m 21 -H 1000000000 -t 10 Friedmanniomyces_simplex_CCFEE_5184.genomic.fa Friedmanniomyces_simplex_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simp_contam
kat comp -m 21 -H 1000000000 -t 10 -o Heterozygosity_F_simp_assembly_reads 'Friedmanniomyces_simplex_?.fastq' Friedmanniomyces_simplex_CCFEE_5184.genomic.fa
