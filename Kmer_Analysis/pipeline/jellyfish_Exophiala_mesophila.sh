#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out jellyfish_Exomeso.log

module load jellyfish
jellyfish count -C -m 25 -s 1000000000 -t 10 Exophiala_mesophila_1.fastq Exophiala_mesophila_2.fastq -o Exophiala_mesophila_k25_reads.jf

jellyfish count -C -m 21 -s 1000000000 -t 10 Exophiala_mesophila_1.fastq Exophiala_mesophila_2.fastq -o Exophiala_mesophila_k21_reads.jf

jellyfish count -C -m 23 -s 1000000000 -t 10 Exophiala_mesophila_1.fastq Exophiala_mesophila_2.fastq -o Exophiala_mesophila_k23_reads.jf

jellyfish count -C -m 19 -s 1000000000 -t 10 Exophiala_mesophila_1.fastq Exophiala_mesophila_2.fastq -o Exophiala_mesophila_k19_reads.jf
