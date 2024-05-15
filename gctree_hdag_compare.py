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
@click.argument('parsimony_forest')
@click.argument('fivemer_mutabilities')
@click.argument('fivemer_substitution')
@click.argument('dnapars_outfile')
@click.argument('meta_file')
@click.argument('abundances')
@click.argument('root_name')
def main(parsimony_forest, fivemer_mutabilities, fivemer_substitution, dnapars_outfile, meta_file, abundances, root_name):

    with open(meta_file, 'r') as fh:
        chain_split = json.loads(fh.read())["l_offset"]

    parsimony_forest = Path(parsimony_forest)
    with open(parsimony_forest, 'rb') as fh:
        forest = pickle.load(fh)
    dag = forest._forest
    parameters = forest.parameters
    trees = gctree.phylip_parse.parse_outfile(dnapars_outfile, abundances, root_name)
    cforests = [gctree.branching_processes.CollapsedForest([gctree.phylip_parse.disambiguate(tree)]) for tree in trees]
    kwargls = [
        (ll_dagfuncs := gctree.branching_processes._ll_genotype_dagfuncs(*forest.parameters), max, dag, ""),
        (mut_funcs := gctree.branching_processes._mutability_dagfuncs(
            splits=[int(chain_split)],
            mutability_file=fivemer_mutabilities,
            substitution_file=fivemer_substitution,
        ), min, dag, ""),
        (placeholder_funcs := hdag.utils.AddFuncDict(
            {
                "start_func": lambda n: 0,
                "edge_weight_func": lambda n1, n2: 0,
                "accum_func": sum,
            },
            name="Whole DAG",
        ), min, dag, ""),
    ]

    ranking_funcs = ll_dagfuncs + mut_funcs

    dag_trees = dag.count_histories()
    best_ll = dag.trim_optimal_weight(**ll_dagfuncs, optimal_func=max)
    best_mutability = dag.trim_optimal_weight(**mut_funcs, optimal_func=min)


    def sign_switch(tup):
        return -tup[0], tup[1]

    best_nondag_tree = min(cforests, key=lambda f: sign_switch(f._forest.optimal_weight_annotate(**ranking_funcs)))

    nondag_best_ll, nondag_best_mutability = best_nondag_tree._forest.optimal_weight_annotate(**ranking_funcs)

    nondag_trees = len(cforests)

    print(f"{parsimony_forest},False,{nondag_trees},{nondag_best_ll},{nondag_best_mutability}")
    print(f"{parsimony_forest},True,{dag_trees},{best_ll},{best_mutability}")

    # with open(output_dir / 'rf_data_counters.p', 'wb') as fh:
    #     fh.write(pickle.dumps(rf_data))

    # print("==========================")
    # print(rf_data)
        

if __name__ == "__main__":
    main()
