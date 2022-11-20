#!/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out jellyfish_Racant.log

module load jellyfish
SAMPFILE=F_endo.txt

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

hostname
date
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read READ
do
  echo "STRAIN is $STRAIN"
  jellyfish count -C -m 25 -s 1000000000 -t 10 ${READ}_1.fastq ${READ}_2.fastq -o {READ}_k25_reads.jf
  jellyfish count -C -m 21 -s 1000000000 -t 10 ${READ}_1.fastq ${READ}_2.fastq -o {READ}_k21_reads.jf
  jellyfish count -C -m 23 -s 1000000000 -t 10 ${READ}_1.fastq ${READ}_2.fastq -o {READ}_k23_reads.jf
  jellyfish count -C -m 19 -s 1000000000 -t 10 ${READ}_1.fastq ${READ}_2.fastq -o {READ}_k19_reads.jf
done
