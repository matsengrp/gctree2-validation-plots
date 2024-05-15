set -eu
# # Uncomment for cluster or server (headless)
# ml libGLU/9.0.2-GCCcore-11.2.0
# ml GLib/2.69.1-GCCcore-11.2.0

TMPDIR=/tmp
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
# eval “$(conda shell.bash hook)”

# Provide to this script the partis output directory containing simulations,
# and an output directory in which to place the gctree outputs for each
# simulation, and an S5F mutability and substitution file
input_directory=$1
output_directory=$2
fivemer_mutabilities=$3
fivemer_substitution=$4

/home/wdumm/miniconda3/envs/partis/bin/python2 /home/wdumm/partis/extract_gctree_inputs.py $input_directory $output_directory

for gctreedir in $output_directory/*; do
    alignment=$gctreedir/all_seqs.fasta
    inference_fasta=$gctreedir/inference.fasta
    output_directory=$gctreedir/gctree
    true_tree=$gctreedir/true_tree.nwk
    chain_split=$(cat $gctreedir/chain_split_igh_len.txt)
    frame1=$(cat $gctreedir/frame1_igh.txt)
    frame2=$(cat $gctreedir/frame2_igk.txt)

    mkdir -p $output_directory
    idmap=$output_directory/idmap.txt
    abundances=$output_directory/abundances.csv
    phylipfile=$(realpath $output_directory/deduplicated.phylip)
    dnaparscfg=$(realpath $output_directory/dnapars.cfg)
    dnaparslog=$(realpath $output_directory/dnapars.log)


    deduplicate $inference_fasta --root naive --idmapfile $idmap --abundance_file $abundances > $phylipfile
    mkconfig $phylipfile dnapars > $dnaparscfg
    currdir=$(pwd)
    # so that dnapars outputs will go where they need to:
    cd $output_directory
    dnapars < $dnaparscfg > $dnaparslog
    cd $currdir

    base_inference_name=$output_directory/gctree_base
    parsimony_forest=$base_inference_name.inference.parsimony_forest.p

    gctree infer $output_directory/outfile $abundances \
        --outbase $base_inference_name \
        --root naive \
        --idmapfile $idmap \
        --chain_split $chain_split \
        --frame $frame1 \
        --frame2 $frame2 \
        --verbose \
        | tee $base_inference_name.log

    mut_inference_name=$output_directory/gctree_mle_mut
    gctree infer $parsimony_forest \
        --outbase $mut_inference_name \
        --root naive \
        --idmapfile $idmap \
        --chain_split $chain_split \
        --frame $frame1 \
        --frame2 $frame2 \
        --mutability $fivemer_mutabilities \
        --substitution $fivemer_substitution \
        --verbose \
        | tee $mut_inference_name.log

    python gctree_benchmark.py $output_directory $parsimony_forest $fivemer_mutabilities $fivemer_substitution $true_tree $alignment $base_inference_name $mut_inference_name $chain_split
done
