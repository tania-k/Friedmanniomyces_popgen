#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 64G --out FastTreeMP.%A.log --time 48:00:00 -p stajichlab

module load fasttree/2.1.11
module unload perl
module unload python
module unload miniconda2
module load miniconda3


IN=Friedmannomyces_Life_Paper.15_taxa.fungi_odb10.aa.fasaln
OUT=$(basename $IN .fasaln).FT.tre
OUTLONG=$(basename $IN .fasaln).FT_long.tre

FastTreeMP -noml < $IN > $OUT
perl PHYling_unified/util/rename_tree_nodes.pl $OUT prefix.tab > $OUTLONG
