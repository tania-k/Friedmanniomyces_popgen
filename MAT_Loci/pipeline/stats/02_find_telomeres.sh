#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 4 --mem 24gb  --out logs/find_telomeres.log

module load parallel

mkdir -p reports
ls cds/*.fa | parallel -j 4 python  scripts_Hiltunen/find_telomeres.py {} \> reports/{/.}.telomere_report.txt
