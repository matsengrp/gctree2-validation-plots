#!/bin/zsh
set -eu
# Modify these prefixes to point to simulations and inference results from the
# simulation pipeline provided with the paper:
sim_prefix=/fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/v6
inference_prefix=/fh/fast/matsen_e/dralph/partis/paired-loci/gct-valid/v6
# If a cluster is not available, change value to 0:
use_cluster=1

conda activate gctree2-validation-plots

wget -O HS5F_Mutability.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Mutability.csv
wget -O HS5F_Substitution.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Substitution.csv

mkdir -p workdir
true_filepaths=workdir/truefilepathsmap.txt
true_trees=workdir/truetreesmap.txt
rm -f $true_filepaths
touch $true_filepaths

# #### Build list of true tree newicks
rm -f $true_trees
touch $true_trees

for gctreepath in $sim_prefix/seed-*/obs-times-(15|20|30|40|50)/simu/selection/simu/event-*/; do
    echo $gctreepath $(./get_true_tree_newick.py $gctreepath/simu_lineage_tree.p) >> $true_trees
done
# ## End building Newicks of true trees

echo Done putting new trees in $true_trees


# Do evaluation:
count=0
for gctreepath in $inference_prefix/seed-*/obs-times-(15|20|30|40|50)/partis/gctree/iclust-*/; do
    # if a cluster is available:
    if [ "$use_cluster" -eq 1 ]; then
        sbatch -c 1 -J bn$count -o cluster_benchmark.log \
            --wrap "\
                python gctree_benchmark_direct.py testdir/$count.p $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/input-seqs.fa $gctreepath/meta.yaml $true_trees $gctreepath/outfile $gctreepath/abundances.csv XnaiveX"
    else
        python gctree_benchmark_direct.py workdir/$count.p $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/input-seqs.fa $gctreepath/meta.yaml $true_trees $gctreepath/outfile $gctreepath/abundances.csv XnaiveX
    fi

    echo $gctreepath $count >> $true_filepaths
    echo $gctreepath $count
    count=$(expr $count + 1)
done


if [ "$use_cluster" -eq 1 ]; then
    # Wait for all cluster jobs to finish (Comment out if no cluster!)
    while [ "$(squeue -u $USER -l | wc -l)" -ne 2 ]; do
      # Sleep for a short interval before checking again
      sleep 30
    done
fi


# Do plotting:
python plot_criterion_comparison.py

echo "Building example tree scatters for all trees in example_sims.txt. To try more from the table above, you can look up their filepaths in " $true_filepaths " then rerun tree_scatter.py"


rm -rf example_sim_scatters
mkdir example_sim_scatters

cat example_sims.txt | while read -r sim; do
    # For each simulation number, find the corresponding line in truefilepathsmap.txt
    gctreepath=$(grep " $sim$" $true_filepaths | awk '{print $1}' | head -n 1)
    
    # Check if a match was found
    if [ -n "$gctreepath" ]; then
        echo "Simulation $sim: $gctreepath"
        python gctree_benchmark_direct.py ignorethisfile.p $gctreepath/gctree.out.inference.parsimony_forest.p HS5F_Mutability.csv HS5F_Substitution.csv $gctreepath/input-seqs.fa $gctreepath/meta.yaml $true_trees $gctreepath/outfile $gctreepath/abundances.csv XnaiveX -a workdir/all_dagtrees_example.p

        rm ignorethisfile.p
        python tree_scatter.py
        mv tree_scatter.pdf example_sim_scatters/$sim.pdf

    else
        echo "Simulation $sim: No matching file path found"
    fi
done
