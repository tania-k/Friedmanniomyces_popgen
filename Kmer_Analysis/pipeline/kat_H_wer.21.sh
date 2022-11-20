#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out logs/H_wer.kat_k21.log

conda activate kat
kat hist -m 21 -H 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k21_kat.jf
kat gcp -m 21 -H 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k21.gcp
kat comp -m 21 -H 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_k21.comp
kat comp -H 1000000000 -t 10 Hortaea_werneckii_1.fastq Hortaea_werneckii_2.fastq -o Hortaea_werneckii_all.comp
