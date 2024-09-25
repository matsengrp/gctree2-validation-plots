#!/home/wdumm/gctree-benchmark-new/gctreeenvpy3.10/bin/python
# (before switching to python 3, the necessary path was
# /home/wdumm/miniconda3/envs/partis/bin/python)
import pickle

import sys

input_pickle = sys.argv[1]

with open(input_pickle, "rb") as fh:
    t = pickle.load(fh)

print(t.write(features=["nuc_seq"], format_root_node=True))
