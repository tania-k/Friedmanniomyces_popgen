#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 64G --out FastTreeMP.%A.log --time 98:00:00 -p intel
CPU=24
module load IQ-TREE/2.1.1
IN=Friedmannomyces_Life_Paper.15_taxa.fungi_odb10.aa.fasaln
PRE=$(basename $IN _taxa.JGI_1086.fasaln)
iqtree -alrt 1000 -bb 1000 -s $IN -pre $PRE -nt AUTO
