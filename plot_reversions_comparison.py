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

simu_reversions = [dc[1]["SimuReversions"] for dc in datacounters]
infer_reversions = [dc[1]["LikelihoodThenContextReversions"] for dc in datacounters]


# Create a figure and axis objects
fig, ax = plt.subplots()


# Plot normalized histogram for real data with bins one unit wide
counts, bins, patches = ax.hist(
    simu_reversions, alpha=0.5, label="Simulated Trees", density=True
)

# # Annotate the bars with the out-degree number and remove x-axis ticks
# ax.set_xticks([])
# for count, bin in zip(counts, bins[:-1]):
#     if count > 0 and int(bin) % 2 == 0:  # Only annotate bars with data
#         ax.text(bin + 0.5, -.002, f'{int(bin)}', ha='center', va='bottom', color='black', fontsize=5, zorder=1000)

# Plot outlined histogram for simulated data with bins one unit wide
ax.hist(
    infer_reversions,
    alpha=0.5,
    label="Inferred Trees",
    density=True,
    histtype="step",
    edgecolor="black",
    bins=bins,
)


# Add labels and title
ax.set_xlabel("Reversion Counts")
ax.set_ylabel("Normalized Frequency")
ax.set_title("Reversion Count Comparison")
# ax.xaxis.labelpad = 20
ax.legend()

# Show plot
fig.savefig("reversion_comparison.pdf")
