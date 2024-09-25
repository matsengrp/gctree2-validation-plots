import polars as pl
import statistics as stat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import numpy as np
from pathlib import Path
import pickle
from counterstats import countermean, countermedian
from matplotlib import colormaps
from collections import defaultdict, Counter

read_data_from = Path("testdir")


counter_summarizers = [
    ("min", lambda c: min(c.keys())),
    ("max", lambda c: max(c.keys())),
    ("mean", countermean),
    ("median", countermedian),
    ("treecount", lambda c: str(sum(c.values()))),
]

dnapars_summarizers = [
    ("min", min),
    ("max", max),
    ("mean", stat.mean),
    ("median", stat.median),
]


def load_counters(cpath):
    with open(cpath, "rb") as fh:
        return pickle.load(fh)


datacounters = [
    (int(p.stem), load_counters(p))
    for p in read_data_from.glob("*")
    if p.stem not in ("truefilepathsmap", "truetreesmap", "all_dagtrees_example")
]


datadict = defaultdict(list)


def simulation_has_context(sim_log_path):
    with open(sim_log_path, "r") as fh:
        for line in fh:
            if "no_context" in line:
                return False
        return True


rf_handlers = [
    (
        lambda raw_rf, infercount, truecount, numleaves: raw_rf
        / (infercount + truecount - 2 * numleaves),
        "Normalized_rf",
    ),
    (lambda raw_rf, infercount, truecount, numleaves: raw_rf, "Raw_rf"),
]

lex_val_names = {
    "best_Likelihood_then_Context": ("BPLikelihoodLogLoss", "ContextLikelihoodLogLoss"),
    "best_Reversions_then_Default": (
        "NaiveReversions",
        "BPLikelihoodLogLoss",
        "ContextLikelihoodLogLoss",
    ),
    "best_Reversions_then_Context": ("NaiveReversions", "ContextLikelihoodLogLoss"),
}

for sim_id, simdict in datacounters:
    counterdict = simdict["WholeDAG"]

    datadict["SimID"].append(sim_id)
    datadict["num_leaves"].append(simdict["NumLeaves"])
    datadict["TrueTreeNumNodes"].append(simdict["TrueTreeNumNodes"])

    num_true_nodes = simdict["TrueTreeNumNodes"]
    num_leaves = simdict["NumLeaves"]

    def process_rf_counter(rf_nodecount_counter, process_func):
        # if this is changed, must also change how the nondag values are
        # computed below...
        newcounter = Counter()
        for (rf, nodecount), count in rf_nodecount_counter.items():
            newcounter.update(
                {process_func(rf, nodecount, num_true_nodes, num_leaves): count}
            )
        return newcounter

    for handler, handler_name in rf_handlers:
        for dtype_name, datacounter in counterdict.items():
            norm_rfcounter = process_rf_counter(datacounter, handler)
            for summary_name, summarizer in counter_summarizers:
                datadict[summary_name + dtype_name + handler_name].append(
                    summarizer(norm_rfcounter)
                )

    dnapars_datafields = simdict["dnaparsTrees"][0]
    dnapars_data = simdict["dnaparsTrees"][1:]
    datadict["dnapars_NumTrees"].append(len(dnapars_data))

    new_dnapars_datafields = (
        tuple([it[1] for it in rf_handlers]) + dnapars_datafields[1:]
    )

    def process_dnapars_tuples(tup):
        rf, infercount = tup[0]
        return (
            tuple(
                handler(rf, infercount, num_true_nodes, num_leaves)
                for handler, _ in rf_handlers
            )
            + tup[1:]
        )

    new_dnapars_data = [process_dnapars_tuples(tup) for tup in dnapars_data]

    for idx, name in enumerate(new_dnapars_datafields):
        for sum_name, summarizer in dnapars_summarizers:
            datadict["dnapars_" + sum_name + name].append(
                summarizer([float(t[idx]) for t in new_dnapars_data])
            )

    for idx, name in enumerate(new_dnapars_datafields):
        for opt in (min, max):
            opt_val = opt(tup[idx] for tup in new_dnapars_data)
            datadict["dnapars_" + opt.__name__ + "_" + name + "_mean_norm_rf"].append(
                stat.mean(
                    tup[0]
                    for tup in new_dnapars_data
                    if math.isclose(tup[idx], opt_val)
                )
            )

    datadict["dnapars_BestMutabilityMaxLikelihood"].append(
        min([(-t[2], t[3]) for t in new_dnapars_data])[1]
    )

    for key, val in simdict["WholeDAGTrimVals"].items():
        if isinstance(val, tuple):
            for valname, val in zip(lex_val_names[key], val):
                datadict[valname + "_from_" + key].append(float(val))
        else:
            datadict["trimval_" + key].append(float(val))

    # Summarize some simulation characteristics:
    sim_inference_path = Path(simdict["InferencePath"])
    with open(sim_inference_path.parent / "abundances.csv", "r") as fh:
        abundances = [int(l.split(",")[-1]) for l in fh]
    datadict["max_abundance"].append(float(max(abundances)))
    datadict["fraction_nonsingleton"].append(
        len([it for it in abundances if it > 1]) / len(abundances)
    )
    datadict["context_dependent"].append(
        simulation_has_context(simdict["SimulationPath"] / "../../../simu.log")
    )


unfiltered_df = pl.DataFrame(datadict)

# Filter for informative abundances
abundance_filter = (pl.col("max_abundance") > 2) & (
    pl.col("fraction_nonsingleton") > 0.1
)
# df = unfiltered_df.filter(abundance_filter)
# Filter for context dependence (I don't think any are non-context-dependent)
df = unfiltered_df.filter(pl.col("context_dependent"))

for faceted_df, facet_name in [
    (df.filter(abundance_filter), "informative_abundance"),
    (df.filter(abundance_filter.not_()), "uninformative_abundance"),
]:

    # Creating the box plot
    fig, ax = plt.subplots(figsize=(8, 8))  # Create a figure and an axes
    columns = [
        "meanWhole DAGNormalized_rf",
        "meanLikelihood_then_ContextNormalized_rf",
        "meanReversions_then_DefaultNormalized_rf",
        "meanReversions_then_ContextNormalized_rf",
        "minWhole DAGNormalized_rf",
    ]
    box = ax.boxplot(
        [faceted_df.select((pl.col(col)).alias("this"))["this"] for col in columns],
        patch_artist=True,
        medianprops={"color": "black"},
    )  # Add a box plot to the axes

    # Set colors from the Dark2 colormap
    dark2_colors = colormaps["Dark2"].colors
    for patch, color in zip(box["boxes"], dark2_colors):
        patch.set_facecolor(color)

    # Add horizontal grid lines at every y-axis tick
    ax.yaxis.grid(True, linestyle="-", which="major", color="lightgrey", alpha=0.7)
    # Optional: Customize your plot
    ax.set_title(
        f"mean RF distance for each criterion, {facet_name}, n_sims={len(faceted_df)}"
    )
    ax.set_ylabel("RF distance")
    col_names = [col[4:] for col in columns]
    col_names[0] = "Parsimony Only\n(Whole hDAG)"
    col_names[-1] = "Best MP Tree Found"
    ax.set_xticklabels(col_names)  # Set the labels for each box plot
    ax.tick_params(axis="x", labelrotation=90)
    fig.tight_layout()  # Adjust layout to make room for labels
    fig_name = f"agg_boxplot_{facet_name}.pdf"
    fig.savefig(fig_name)
    print(fig_name)


df = df.filter(abundance_filter)
# Using the object-oriented interface
fig, ax = plt.subplots()


# Add grid lines
ax.grid(True, which="both", linestyle="--", linewidth=0.5)

_x_y_data = ["dnapars_minNormalized_rf", "minWhole DAGNormalized_rf"]

xy = [df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data]
min_val = 0
max_val = max(xy[0].max(), xy[1].max())

scatter = ax.scatter(
    *xy,
    edgecolor="black",  # Black border around each point
)
# Add y=x line
ax.plot(
    [min_val, max_val], [min_val, max_val], "k--", alpha=0.3
)  # 'k--' for black dashed line, alpha for transparency

ax.set_title("Best Tree Improvement due to hDAG")
ax.set_xlabel("Dnapars Best Tree RF Distance")
ax.set_ylabel("hDAG Best Tree RF Distance")
fig_name = "minRF_Scatter.pdf"
fig.savefig(fig_name)
print(fig_name)

# Calculate RF distance improvements and prepare data for box plot
improvements = {"BPLikelihoodLogLoss": [], "ContextLikelihoodLogLoss": []}
for opt, criterion in [
    ("min", "BPLikelihoodLogLoss"),
    ("min", "ContextLikelihoodLogLoss"),
]:
    _x_y_data = [
        f"dnapars_{opt}_{criterion}_mean_norm_rf",
        f"mean{criterion}Normalized_rf",
    ]
    x, y = [df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data]
    improvement = y - x
    improvements[criterion] = improvement

improvements["max improvement"] = (
    df["minWhole DAGNormalized_rf"] - df["dnapars_minNormalized_rf"]
)

# Creating the box plot
fig, ax = plt.subplots()

# Adding grid lines
ax.grid(True, which="both", linestyle="--", linewidth=0.5)

# Create box plot
boxplot = ax.boxplot(
    improvements.values(), labels=improvements.keys(), patch_artist=True
)

# Add titles and labels
ax.set_title("RF Distance Improvement in hDAG Relative to Dnapars")
ax.set_ylabel("RF Distance Improvement")

# Optionally, add colors or other styles to the box plot
colors = colormaps["Dark2"].colors
for patch, color in zip(boxplot["boxes"], colors):
    patch.set_facecolor(color)

fig_name = "RF_distance_improvement_boxplot.pdf"
fig.savefig(fig_name)
print(fig_name)

for criterion in ["BPLikelihoodLogLoss", "ContextLikelihoodLogLoss"]:
    # Using the object-oriented interface
    fig, ax = plt.subplots()

    # Add grid lines
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    _x_y_data = [
        f"dnapars_max_{criterion}_mean_norm_rf",
        f"mean{criterion}Normalized_rf",
    ]
    xy = [df.select((pl.col(col)).alias("this"))["this"] for col in _x_y_data]
    min_val = 0
    max_val = max(xy[0].max(), xy[1].max())

    scatter = ax.scatter(
        *xy,
        edgecolor="black",  # Black border around each point
    )
    # Add y=x line
    ax.plot(
        [min_val, max_val], [min_val, max_val], "k--", alpha=0.3
    )  # 'k--' for black dashed line, alpha for transparency

    ax.set_title(f"mean RF improvement with hDAG, trimmed by {criterion}")
    ax.set_xlabel("Dnapars Best Tree RF Distance")
    ax.set_ylabel("hDAG Best Tree RF Distance")
    fig_name = (
        f"{criterion.replace(' ', '_').replace('.','')}_hdag_improve_RF_Scatter.pdf"
    )
    fig.savefig(fig_name)
    print(fig_name)

# df = pl.read_csv('gctree_hdag_compare_output.csv', dtypes={"NumTreesRanked": pl.String})
# non_hdag_df = df.filter(pl.col('FromHistoryDAG') == False)
# hdag_df = df.filter(pl.col('FromHistoryDAG'))
# hdag_df = hdag_df.rename(lambda name: name if name == "ParsimonyForestPath" else "hDAG_" + name)
# combined_df = hdag_df.join(non_hdag_df, on="ParsimonyForestPath")
#
hdag_comparison_df = pl.DataFrame(
    [
        pl.Series(
            name="Log_hDAG_NumTreesRanked",
            values=[math.log(int(x)) for x in df["treecountWhole DAGRaw_rf"]],
        ),
        pl.Series(
            name="Log_dnapars_NumTrees",
            values=[math.log(int(x)) for x in df["dnapars_NumTrees"]],
        ),
        pl.Series(
            name="RF_Improvement",
            values=(df["minWhole DAGRaw_rf"] < df["dnapars_minRaw_rf"]),
        ),
        df["trimval_best_ContextLikelihoodLogLoss"],
        df["dnapars_minContextLikelihoodLogLoss"],
        df["trimval_best_BPLikelihoodLogLoss"],
        df["dnapars_minBPLikelihoodLogLoss"],
        df["SimID"],
    ]
)

hdag_comparison_df = hdag_comparison_df.with_columns(
    (
        -pl.col("trimval_best_ContextLikelihoodLogLoss")
        + pl.col("dnapars_minContextLikelihoodLogLoss")
    ).alias("Context Likelihood Improvement"),
    (
        -pl.col("trimval_best_BPLikelihoodLogLoss")
        + pl.col("dnapars_minBPLikelihoodLogLoss")
    ).alias("Branching Process Likelihood Improvement"),
    (pl.col("Log_hDAG_NumTreesRanked") - pl.col("Log_dnapars_NumTrees")).alias(
        "Trees Ranked Improvement"
    ),
)


# Create a plot using the object-oriented interface
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)


ax1.yaxis.grid(
    True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
)
ax2.yaxis.grid(
    True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
)
ax1.xaxis.grid(
    True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
)
ax2.xaxis.grid(
    True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
)
ax1.axhline(
    y=1,
    color="lightgrey",
    linestyle="-",
    linewidth=1.5,
    label="Reference (y=1)",
    zorder=3,
)
ax2.axhline(
    y=1,
    color="lightgrey",
    linestyle="-",
    linewidth=1.5,
    label="Reference (y=1)",
    zorder=3,
)
ax1.axvline(
    x=1,
    color="lightgrey",
    linestyle="-",
    linewidth=1.5,
    label="Reference (x=1)",
    zorder=3,
)
ax2.axvline(
    x=1,
    color="lightgrey",
    linestyle="-",
    linewidth=1.5,
    label="Reference (x=1)",
    zorder=3,
)
# Plotting the data
ax1.scatter(
    hdag_comparison_df["Trees Ranked Improvement"].exp(),
    hdag_comparison_df["Branching Process Likelihood Improvement"].exp(),
    color=colormaps["Dark2"].colors[0],
    alpha=0.5,
    s=20,
    linewidths=0,
    label="Branching Process Likelihood Improvement",
    zorder=4,
)
ax2.scatter(
    hdag_comparison_df["Trees Ranked Improvement"].exp(),
    hdag_comparison_df["Context Likelihood Improvement"].exp(),
    facecolors=colormaps["Dark2"].colors[1],
    alpha=0.5,
    s=20,
    linewidths=0,
    label="Context Likelihood Improvement",
    zorder=4,
)

# Setting log scale for axes
ax1.set_xscale("log")
ax2.set_yscale("log")

# Adding labels and title
ax2.set_xlabel("Fold Increase in Trees Ranked")
ax1.set_ylabel("Branching Process\nLikelihood Improvement")
ax2.set_ylabel("Context Likelihood\nImprovement")

# Adding legend
filename = "hdag_comparison.pdf"
fig.savefig(filename)
print(filename)

# Create a plot using the object-oriented interface
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# Apply grid settings and reference lines
for ax in (ax1, ax2):
    ax.yaxis.grid(
        True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
    )
    ax.xaxis.grid(
        True, linestyle="-", which="major", color="lightgrey", alpha=0.7, zorder=3
    )
    ax.axhline(y=1, color="lightgrey", linestyle="-", linewidth=1.5, zorder=3)
    ax.axvline(x=1, color="lightgrey", linestyle="-", linewidth=1.5, zorder=3)


# Helper function to categorize data
def categorize_data(df, y_col):
    df = df.with_columns((df[y_col] > 0).alias("category"))
    proportion_above = df.filter(pl.col("category")).height / df.height
    return df, proportion_above


# Categorize data
hdag_comparison_df, prop_above1 = categorize_data(
    hdag_comparison_df, "Branching Process Likelihood Improvement"
)
hdag_comparison_df2, prop_above2 = categorize_data(
    hdag_comparison_df, "Context Likelihood Improvement"
)


# Function to plot categorized data
def plot_categorized_data(ax, df, x_col, y_col, color_base, label_base):
    for category, color in zip(
        [True, False], [color_base, mcolors.to_rgba(color_base, alpha=0.7)]
    ):
        subset = df.filter(pl.col("category") == category)
        linewidths = 1 if category else 0
        ax.scatter(
            subset[x_col].exp().to_numpy(),
            subset[y_col].exp().to_numpy(),
            color=color,
            alpha=0.5,
            s=20,
            edgecolor="#000000",
            linewidths=linewidths,
            label=f"{'' if category else 'No '} Improvement: {subset.height / df.height:.2%}",
            zorder=4,
        )


# Plotting the data
plot_categorized_data(
    ax1,
    hdag_comparison_df,
    "Trees Ranked Improvement",
    "Branching Process Likelihood Improvement",
    colormaps["Dark2"].colors[3],
    "Branching Process Likelihood Improvement",
)
plot_categorized_data(
    ax2,
    hdag_comparison_df2,
    "Trees Ranked Improvement",
    "Context Likelihood Improvement",
    colormaps["Dark2"].colors[3],
    "Context Likelihood Improvement",
)

# Adding labels and title
ax2.set_xlabel("Fold Increase in Trees Ranked")
ax1.set_ylabel("Branching Process\nLikelihood Improvement")
ax2.set_ylabel("Context Likelihood\nImprovement")

# Setting log scale for axes
ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")


# Adding legend
ax1.legend()
ax2.legend()

# Save the figure
filename = "hdag_comparison_faceted.pdf"
fig.savefig(filename)
print(filename)


fig, axs = plt.subplots(3, 1, figsize=(8, 12))

# Histogram for Context Likelihood Improvement
axs[0].hist(
    hdag_comparison_df["Context Likelihood Improvement"].exp(),
    bins=10,
    color="blue",
    alpha=0.7,
)
axs[0].set_title("Context Likelihood Improvement")
axs[0].set_xlabel("Improvement")
axs[0].set_ylabel("Frequency")

# Histogram for Branching Process Likelihood Improvement
axs[1].hist(
    hdag_comparison_df["Branching Process Likelihood Improvement"].exp(),
    bins=10,
    color="red",
    alpha=0.7,
)
axs[1].set_title("Branching Process Likelihood Improvement")
axs[1].set_xlabel("Improvement")
axs[1].set_ylabel("Frequency")

# Histogram for Trees Ranked Improvement
axs[2].hist(
    hdag_comparison_df["Trees Ranked Improvement"].exp(),
    bins=10,
    color="green",
    alpha=0.7,
)
axs[2].set_title("Trees Ranked Improvement")
axs[2].set_xlabel("Improvement")
axs[2].set_ylabel("Frequency")

plt.tight_layout()
filename = "hdag_comparison_histograms.pdf"
fig.savefig(filename)
print(filename)


# fig, axs = plt.subplots(3, 1, figsize=(8, 12))
#
# hdag_comparison_df = hdag_comparison_df.sort("Trees Ranked Improvement")
# # Plotting rank plots for each improvement metric
# for i, key in enumerate(["Context Likelihood Improvement", "Branching Process Likelihood Improvement", "Trees Ranked Improvement"]):
#     # Sorting the data and obtaining ranks
#     # sorted_data = np.sort(hdag_comparison_df[key].exp())
#     sorted_data = hdag_comparison_df[key].exp()
#     ranks = np.arange(1, len(sorted_data) + 1)
#
#     # Plotting
#     axs[i].scatter(ranks, sorted_data, marker='o', linestyle='-', edgecolors='k')
#     axs[i].set_title(f'{key} (Ranked)')
#     axs[i].set_xlabel('Rank')
#     axs[i].set_ylabel('Improvement')
#     axs[i].set_yscale('log')
#     axs[i].axhline(y=1, color='lightgrey', linestyle='-', linewidth=1.5, label='Reference (y=1)')
#     # axs[i].set_ylim(bottom=0)
#     # y_ticks = np.unique(np.append([1], axs[i].get_yticks()))
#     # axs[i].set_yticks(y_ticks)
#     axs[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7)
#
# plt.tight_layout()
# filename = "hdag_comparison_rankplots.pdf"
# fig.savefig(filename)
# print(filename)
# # Add a color bar as a legend for the colors
# cbar = fig.colorbar(scatter, ax=ax)
# cbar.set_label('Ratio of Best Mutability')

# fig.savefig("hdag_comparison.pdf")


# # Using the object-oriented interface
# fig, ax = plt.subplots()

# scatter = ax.scatter(
#     combined_df["RatioBestMutability"],
#     combined_df["LogRatioBestLikelihood"],
#     c=combined_df["Log_NumTreesRanked"],
#     edgecolor='black',  # Black border around each point
# )

# ax.set_xlabel('Ratio of Best Mutability')
# ax.set_ylabel('Log Ratio of Best Likelihood')

# # Add a color bar as a legend for the colors
# cbar = fig.colorbar(scatter, ax=ax)
# cbar.set_label('Log Ratio Trees Ranked')

# fig.savefig("hdag_comparison_mut_v_likelihood.pdf")

# # Assuming combined_df is already defined and contains the necessary columns
# x_labels = np.arange(len(combined_df))  # Qualitative x-axis based on the number of rows

# # Create a figure with three subplots, sharing the x-axis
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 12))

# # Scatter plot for log number of trees ranked
# ax1.scatter(x_labels, combined_df["LogRatioTreesRanked"])
# ax1.axhline(y=0, color='gray', linestyle='--')  # Reference line at y=0
# ax1.set_ylabel('Log Number of Trees Ranked')

# # Scatter plot for log ratio of best likelihood
# ax2.scatter(x_labels, combined_df["LogRatioBestLikelihood"])
# ax2.axhline(y=0, color='gray', linestyle='--')  # Reference line at y=0
# ax2.set_ylabel('Log Ratio of Best Likelihood')

# # Scatter plot for ratio of best mutability
# ax3.scatter(x_labels, combined_df["RatioBestMutability"])
# ax3.axhline(y=1, color='gray', linestyle='--')  # Reference line at y=1
# ax3.set_ylabel('Ratio of Best Mutability')

# fig.tight_layout()
# fig.savefig("stacked_scatter_plots.pdf")


# ### For choosing a tree_scatter replicate:

tree_scatter_test_df = (
    hdag_comparison_df.sort(
        "Trees Ranked Improvement", "Log_dnapars_NumTrees", "Log_hDAG_NumTreesRanked"
    )
    .filter(pl.col("Log_dnapars_NumTrees") > 2)
    .filter(pl.col("Context Likelihood Improvement") > 0)[
        [
            "SimID",
            "Log_hDAG_NumTreesRanked",
            "Log_dnapars_NumTrees",
            "RF_Improvement",
            "Context Likelihood Improvement",
            "Branching Process Likelihood Improvement",
        ]
    ]
    .sort("RF_Improvement")
)

print("Choose one of these replicates to run tree_scatter.py on:")
pl.Config.set_tbl_rows(100)
print(tree_scatter_test_df.tail(100))

print("All those with RF distance improvement will be written to example_paths.txt.")
simids = list(tree_scatter_test_df.filter(pl.col("RF_Improvement"))["SimID"])

with open("example_sims.txt", "w") as fh:
    for name in simids:
        print(name, file=fh)

# Sim 290 was used in the paper
