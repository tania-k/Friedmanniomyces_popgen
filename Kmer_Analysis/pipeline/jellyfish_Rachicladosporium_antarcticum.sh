#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out jellyfish_Racant.log

module load jellyfish
jellyfish count -C -m 25 -s 1000000000 -t 10 Rachicladosporium_antarcticum_1.fastq Rachicladosporium_antarcticum_2.fastq -o Rachicladosporium_antarcticum_k25_reads.jf

jellyfish count -C -m 21 -s 1000000000 -t 10 Rachicladosporium_antarcticum_1.fastq Rachicladosporium_antarcticum_2.fastq -o Rachicladosporium_antarcticum_k21_reads.jf

jellyfish count -C -m 23 -s 1000000000 -t 10 Rachicladosporium_antarcticum_1.fastq Rachicladosporium_antarcticum_2.fastq -o Rachicladosporium_antarcticum_k23_reads.jf

jellyfish count -C -m 19 -s 1000000000 -t 10 Rachicladosporium_antarcticum_1.fastq Rachicladosporium_antarcticum_2.fastq -o Rachicladosporium_antarcticum_k19_reads.jf
