This directory contains scripts for benchmarking the newest version of gctree,
using Duncan's bcr-phylo simulations.

# Environment Setup:
Running the scripts in this repo requires a python environment with `gctree`
and other dependencies installed. The dependencies are listed in
`requirements.txt` and are all installable with pip. However, the conda
environment file `environment.yml` is also provided.

```
conda env create -f environment.yml
conda activate gctree2-validation-plots
```


# Producing plots:

To produce the plots, you must first edit `gctree_benchmark_direct_run.sh` to
point it toward the simulated data and inference outputs. This may be done by
editing the variables `sim_prefix` and `inference_prefix` at the top of the
script.
The values for `sim_prefix` and `inference_prefix` specify where the outputs
from running the simulation and inference code, which is provided with the
paper.

To produce the plots, run `gctree_benchmark_direct_run.sh` in the environment
set up as described above.
This script will send jobs to a cluster using `sbatch`. If you do not have
access to a cluster using `sbatch`, the `use_cluster` variable at the top of
the script will need to be set to 0.
The script will not finish until all sbatch jobs for your user are finished
running, so be sure all unrelated cluster jobs are finished before running this
script.

# Plot outputs:

Figure 4 from the paper is the file `hdag_comparison_faceted.pdf`.
A variety of plots similar to Figure 5 will be placed in the directory
`example_sim_scatters`. For the paper, we selected one which showed the
following:
* improvement in RF distance in history sDAG trees compared to dnapars trees,
* many more history sDAG trees than dnapars trees, and
* qualitative correlation between RF distance improvement and the two
    likelihoods.
