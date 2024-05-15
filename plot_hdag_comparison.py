import polars as pl
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

df = pl.read_csv('gctree_hdag_compare_output.csv', dtypes={"NumTreesRanked": pl.String})
non_hdag_df = df.filter(pl.col('FromHistoryDAG') == False)
hdag_df = df.filter(pl.col('FromHistoryDAG'))
hdag_df = hdag_df.rename(lambda name: name if name == "ParsimonyForestPath" else "hDAG_" + name)
combined_df = hdag_df.join(non_hdag_df, on="ParsimonyForestPath")

combined_df = combined_df.with_columns(pl.Series(name="Log_hDAG_NumTreesRanked", values=[math.log(int(x)) for x in combined_df["hDAG_NumTreesRanked"]]))
combined_df = combined_df.with_columns(pl.Series(name="Log_NumTreesRanked", values=[math.log(int(x)) for x in combined_df["NumTreesRanked"]]))

combined_df = combined_df.with_columns(
    (pl.col("Log_hDAG_NumTreesRanked") - pl.col("Log_NumTreesRanked")).alias("LogRatioTreesRanked"),
    (pl.col("hDAG_BestLikelihood") - pl.col("BestLikelihood")).alias("LogRatioBestLikelihood"),
    (pl.col("BestMutability") / pl.col("hDAG_BestMutability")).alias("RatioBestMutability"),# order swapped so that >1 is improvement with hdag
)

ratio_min, ratio_max = combined_df["RatioBestMutability"].min(), combined_df["RatioBestMutability"].max()
vmin, vmax = min(ratio_min, 1/ratio_max), max(ratio_max, 1/ratio_min)

combined_df = combined_df.sort(pl.col("LogRatioTreesRanked"))



# # Using the object-oriented interface
# fig, ax = plt.subplots()

# scatter = ax.scatter(
#     combined_df["Log_NumTreesRanked"],
#     combined_df["LogRatioBestLikelihood"],
#     c=combined_df["RatioBestMutability"],
#     cmap='PiYG',  # You can choose a colormap that suits your data
#     edgecolor='black',  # Black border around each point
#     vmin=vmin,
#     vmax=vmax,
# )

# ax.set_title('Improvement in best tree after hDAG')
# ax.set_xlabel('Tree count (log ratio > 0 means improved)')
# ax.set_ylabel('Likelihood (log ratio > 0 means improved)')

# # Add a color bar as a legend for the colors
# cbar = fig.colorbar(scatter, ax=ax)
# cbar.set_label('Mutability Parsimony (> 1 means improved)')

# fig.savefig("hdag_comparison.pdf")


# Using the object-oriented interface
fig, ax = plt.subplots()


# Add grid lines
ax.grid(True, which='both', linestyle='--', linewidth=0.5)

scatter = ax.scatter(
    combined_df["RatioBestMutability"],
    combined_df["LogRatioBestLikelihood"],
    c=combined_df["Log_NumTreesRanked"],
    edgecolor='black',  # Black border around each point
    cmap=LinearSegmentedColormap.from_list("not_viridis", ["white", "purple"]),  # You can choose a colormap that suits your data
)

ax.set_title('Improvement in best tree after hDAG')
ax.set_xlabel('Mutability Parsimony (> 1 means improved)')
ax.set_ylabel('Likelihood (log ratio > 0 means improved)')

# Add a color bar as a legend for the colors
cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label('Tree count (log ratio > 0 means improved)')

fig.savefig("hdag_comparison_mut_v_likelihood.pdf")

# Assuming combined_df is already defined and contains the necessary columns
x_labels = np.arange(len(combined_df))  # Qualitative x-axis based on the number of rows

# Create a figure with three subplots, sharing the x-axis
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 12))

# Scatter plot for log number of trees ranked
ax1.scatter(x_labels, combined_df["LogRatioTreesRanked"])
ax1.axhline(y=0, color='gray', linestyle='--')  # Reference line at y=0
ax1.set_ylabel('Tree count (log ratio > 0 means more trees ranked')

# Scatter plot for log ratio of best likelihood
ax2.scatter(x_labels, combined_df["LogRatioBestLikelihood"])
ax2.axhline(y=0, color='gray', linestyle='--')  # Reference line at y=0
ax2.set_ylabel('Likelihood (log ratio > 0 means improved)')

# Scatter plot for ratio of best mutability
ax3.scatter(x_labels, combined_df["RatioBestMutability"])
ax3.axhline(y=1, color='gray', linestyle='--')  # Reference line at y=1
ax3.set_ylabel('Mutability Parsimony (ratio > 1 means improved)')

fig.suptitle('Improvement in best tree after hDAG, by simulation')
fig.tight_layout()
fig.savefig("stacked_scatter_plots.pdf")
