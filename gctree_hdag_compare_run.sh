#!/bin/bash
set -eu

source /home/wdumm/gctree-benchmark-new/gctreeenvpy3.10/bin/activate
outfile=gctree_hdag_compare_output.csv
echo "ParsimonyForestPath,FromHistoryDAG,NumTreesRanked,BestLikelihood,BestMutability" > $outfile

for gctreepath in /fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/*/obs-times-*/*/partis/gctree/iclust-*/; do
    python gctree_hdag_compare.py $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/outfile $gctreepath/meta.yaml $gctreepath/abundances.csv XnaiveX >> $outfile
    echo $gctreepath

done
