#!/bin/zsh
set -eu

source /home/wdumm/gctree-benchmark-new/gctreeenvpy3.10/bin/activate

mkdir -p testdir
rm -f testdir/*.p
true_filepaths=testdir/truefilepathsmap.txt
true_trees=testdir/truetreesmap.txt
rm -f $true_filepaths
touch $true_filepaths

# #### Build list of true tree newicks
rm -f $true_trees
touch $true_trees

# # This was the old glob pattern:
# for gctreepath in /fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/*/obs-times-*/*/simu/selection/simu/event-*/; do

for gctreepath in /fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/v4/seed-*/obs-times-(15|20|30|40|50)/simu/selection/simu/event-*/; do
    echo $gctreepath $(./get_true_tree_newick.py $gctreepath/simu_lineage_tree.p) >> $true_trees
done
# ## End building Newicks of true trees

echo Done putting new trees in $true_trees

count=0
# for gctreepath in /fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/*/obs-times-*/*/partis/gctree/iclust-*/; do

for gctreepath in /fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/v4/seed-*/obs-times-(15|20|30|40|50)/partis/gctree/iclust-*/; do
    python gctree_benchmark_direct.py testdir/$count.p $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/input-seqs.fa $gctreepath/meta.yaml $true_trees $gctreepath/outfile $gctreepath/abundances.csv XnaiveX

    echo $gctreepath $count >> $true_filepaths
    echo $gctreepath $count
    count=$(expr $count + 1)

done

python plot_criterion_comparison.py

gctreepath=/fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/v4/seed-0/obs-times-40/partis/gctree/iclust-46/
python gctree_benchmark_direct.py ignorethisfile.p $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/input-seqs.fa $gctreepath/meta.yaml $true_trees $gctreepath/outfile $gctreepath/abundances.csv XnaiveX -a testdir/all_dagtrees_example.p

rm ignorethisfile.p
python tree_scatter.py
