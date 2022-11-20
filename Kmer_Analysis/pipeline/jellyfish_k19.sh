#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out jellyfish_k19.log

module load jellyfish
jellyfish count -C -m 19 -s 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_endolithicus_2.fastq -o Friedmanniomyces_endolithicus_k19_reads.jf

jellyfish count -C -m 19 -s 1000000000 -t 10 Friedmanniomyces_endolithicus_1.fastq Friedmanniomyces_simplex_2.fastq -o Friedmanniomyces_simplex_k19_reads.jf
