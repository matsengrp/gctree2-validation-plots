import polars as pl
import statistics as stat
import matplotlib.pyplot as plt
import matplotlib
import math
import numpy as np
from pathlib import Path
import pickle
from counterstats import countermean, countermedian
from matplotlib import colormaps
from collections import defaultdict, Counter

with open('testdir/all_dagtrees_example.p', 'rb') as fh:
    data = pickle.load(fh)

n_leaves, n_nodes_true_tree, dag_data, dnapars_data = data["NumLeaves"], data["TrueTreeNumNodes"], data["dagtrees"], data["dnaparsTrees"]
dag_df = pl.DataFrame({colname: [it[i] for it in dag_data[1:]] for i, colname in enumerate(dag_data[0])})
dnapars_df = pl.DataFrame({colname: [it[i] for it in dnapars_data[1:]] for i, colname in enumerate(dnapars_data[0])})


fig, ax = plt.subplots()



dag_df = dag_df.with_columns(pl.Series("normalized_rf", [rf / (ncount + n_nodes_true_tree - (2 * n_leaves)) for rf, ncount in dag_df["RootedRF_and_nodecount"]]))
dnapars_df = dnapars_df.with_columns(pl.Series("normalized_rf", [rf / (ncount + n_nodes_true_tree - (2 * n_leaves)) for rf, ncount in dnapars_df["RootedRF_and_nodecount"]]))

best_dnapars_rf = dnapars_df["normalized_rf"].min()

# Add grid lines
ax.grid(True, which='both', linestyle='--', linewidth=0.5, zorder=1)

_x_y_data = ["BPLikelihoodLogLoss", "ContextLikelihoodLogLoss"]

scatter_titles = ["History sDAG trees", "Dnapars trees"]

dag_not_better_df = dag_df.filter(pl.col("normalized_rf") >= best_dnapars_rf)
dag_better_df = dag_df.filter(pl.col("normalized_rf") < best_dnapars_rf)

ax.scatter(
    *[-dag_not_better_df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data],
    # edgecolor='black',  # Black border around each point
    linewidths=0,
    color=matplotlib.colors.colorConverter.to_rgba(colormaps["Dark2"].colors[0], alpha=.3),
    # alpha=.6,
    marker='.',
    s=70,
    label="History sDAG trees",
    zorder=2,
)

ax.scatter(
    *[-dnapars_df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data],
    # edgecolor='black',  # Black border around each point
    # facecolors='none',
    linewidths=0,
    color=matplotlib.colors.colorConverter.to_rgba(colormaps["Dark2"].colors[1], alpha=1),
    # alpha=.9,
    marker='.',
    label="Dnapars trees",
    zorder=2,
    s=70,
)

ax.scatter(
    *[-dag_better_df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data],
    edgecolor='black',  # Black border around each point
    linewidths=.7,
    color=matplotlib.colors.colorConverter.to_rgba(colormaps["Dark2"].colors[0], alpha=.3),
    # alpha=.6,
    marker='.',
    s=70,
#     linewidths=0,
#     color=colormaps["Dark2"].colors[1],
#     alpha=.6,
#     marker='.',
    label="History sDAG trees (improved RF-distance)",
    zorder=2,
)


ax.set_title('Dnapars and History sDAG Tree Likelihoods')
ax.set_xlabel('Branching Process Log-Likelihood')
ax.set_ylabel('Poisson Context Log-Likelihood')
ax.legend()
fig_name = "tree_scatter.pdf"
fig.savefig(fig_name)
print(fig_name)
