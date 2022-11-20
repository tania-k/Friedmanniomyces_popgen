#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out compare_sharevsuniqe.log


kat plot spectra-mx -i -r 0 -x 300 -o compare_Friedmanniomyces.png  compare_Friedmanniomyces-main.mx
