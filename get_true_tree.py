"""Computes the ranking criteria for the best ranked tree in the history DAG, as well as for the best ranked tree in only the trees explicitly
found by dnapars. Prints to stdout the following comma-separated columns:
ParsimonyForestPath,FromHistoryDAG,NumTreesRanked,BestLikelihood,BestMutability

Note that trees are ranked lexicographically, maximizing likelihood, then minimizing mutability parsimony. That is, BestMutability is the best mutability
parsimony score of a tree maximizing likelihood, but not necessarily the best mutability parsimony of any tree being ranked.
"""
from collections import Counter
import json
import gctree
import pickle
import historydag as hdag
from pathlib import Path
import ete3
import click

# mutability file
# substitution file
# parsimony forest

@click.command()
@click.argument('input_sequences')
def main(input_sequences):

    input_sequences = Path(input_sequences)
    input_fasta = hdag.utils.load_fasta(input_sequences)
    simu_path = input_sequences.parent / "../../../simu/selection/simu/"
    possible_paths = list(simu_path.glob("event-*"))

    def evaluate_paths(possible_paths):
        for path in possible_paths:
            candidate_fasta = hdag.utils.load_fasta(path / "simu.fasta")
            for key in input_fasta:
                if key != "XnaiveX":
                    if key[0:-15] in candidate_fasta:
                        return path
        
    matched_path = evaluate_paths(possible_paths)
    if matched_path is None:
        raise RuntimeError("Could not match simulation with inference")

    simu_tree = ete3.Tree(newick=(matched_path / 'simu.nwk').as_posix(), format=1)
    simu_alignment = hdag.utils.load_fasta(matched_path / 'simu.fasta')

    # build conversion function which removes padding in simu_alignment
    # sequences to get sequences like the ones in input_fasta


    converted_simu_alignment = {key: val.replace('N', '') for key, val in simu_alignment.items()}


    # Double check that converted sequences are the same for leaves (also need
    # to convert leaf name keys)
    assert {key: val for key, val in converted_simu_alignment.items() if 'leaf' in key} == {key[0:-15]: val for key, val in input_fasta.items() if key != 'XnaiveX'}

    # put sequences from simu_alignment on sim_tree
    for node in simu_tree.traverse():
        node.add_feature("sequence", converted_simu_alignment[node.name])

    # Do any collapsing of simu_tree that may be necessary for comparison with
    # hdag
    seq_counter = Counter((node.sequence for node in simu_tree.iter_leaves()))
    for node in list(simu_tree.traverse()):
        node.add_feature("abundance", seq_counter[node.sequence])

    # and remove duplicate leaves, if possible...
    def delete_duplicates(tree):
        to_delete = []
        visited = set()
        for node in sorted(tree.iter_leaves(), key=lambda n: -n.abundance):
            if node.sequence in visited:
                to_delete.append(node)
            else:
                visited.add(node.sequence)
        num_deleted = len(to_delete)
        for node in to_delete:
            node.delete(prevent_nondicotomic=False)
        return num_deleted

    while delete_duplicates(simu_tree) > 0:
        continue

    # Remove unifurcations:
    to_delete = [n for n in simu_tree.iter_descendants() if len(n.children) == 1]
    for node in to_delete:
        node.delete(prevent_nondicotomic=False)

    cforest = gctree.CollapsedForest([simu_tree])
    return cforest._forest

    
        

if __name__ == "__main__":
    main()
