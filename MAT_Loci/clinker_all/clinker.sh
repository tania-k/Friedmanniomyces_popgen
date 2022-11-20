#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 16gb -J clinker --out clinker.%A.log -p batch

hostname
MEM=64
CPU=$SLURM_CPUS_ON_NODE

source activate clinker

clinker *.gbk -i 0.23 -p Fried_ALL.html
