#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out jellyfish_Hortwer.log

module load jellyfish
jellyfish count -C -m 25 -s 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k25_reads.jf

jellyfish count -C -m 21 -s 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k21_reads.jf

jellyfish count -C -m 23 -s 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k23_reads.jf

jellyfish count -C -m 19 -s 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_reads.jf
