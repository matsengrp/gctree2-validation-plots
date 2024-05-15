from collections import Counter
import gctree
import pickle
import historydag as hdag
from pathlib import Path
import click
import ete3

# mutability file
# substitution file
# parsimony forest

@click.command()
@click.argument('output_dir')
@click.argument('parsimony_forest')
@click.argument('fivemer_mutabilities')
@click.argument('fivemer_substitution')
@click.argument('true_tree')
@click.argument('alignment')
@click.argument('base_inference_name')
@click.argument('mut_inference_name')
@click.argument('chain_split')
def main(output_dir, parsimony_forest, fivemer_mutabilities, fivemer_substitution, true_tree, alignment, base_inference_name, mut_inference_name, chain_split):
    parsimony_forest = Path(parsimony_forest)
    output_dir = Path(output_dir)
    with open(parsimony_forest, 'rb') as fh:
        forest = pickle.load(fh)
    dag = forest._forest
    abundances = {n.label.sequence: n.label.abundance for n in dag.get_leaves()}
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



    fasta = hdag.utils.load_fasta(alignment)
    ##TODO remove line when ambiguities fixed (and in shell script)
    #fasta = {name: seq.replace('N', 'T') for name, seq in fasta.items()}

    tree = ete3.Tree(true_tree, format=1, quoted_node_names=True)
    example_forest_tree = next(iter(forest))
    abundances = {n.sequence: n.abundance for n in example_forest_tree.tree.traverse()}
    for node in tree.traverse():
        node.add_feature("sequence", fasta[node.name])
        if node.is_leaf():
            try:
                node.add_feature("abundance", abundances[node.sequence])
            except KeyError as e:
                print(node.name)
                raise e
        else:
            node.add_feature("abundance", 0)
    # Should we collapse the tree if we use historydag for comparison
    # directly?? (TODO: constructor will fail on trees with convergent
    # evolution resulting in non-unique leaf sequences)
    ctree = gctree.branching_processes.CollapsedTree(tree)

    # test that abundances are correct after collapsing
    ctree_abundances = {n.sequence: n.abundance for n in ctree.tree.traverse()}
    example_abundances = {n.sequence: n.abundance for n in example_forest_tree.tree.traverse()}
    for seq, abundance in ctree_abundances.items():
        if seq in example_abundances:
            assert abundance == example_abundances[seq]
        else:
            assert abundance == 0

    # Gather tree comparison data

    rf_data = {}
    mrca_data = {}
    for kwargs, opt_func, dag, extra_name in kwargls:
        trimmed_dag = dag.copy()
        trimmed_dag.trim_optimal_weight(**kwargs, optimal_func=opt_func)
        trimmed_forest = forest._trimmed_self(trimmed_dag)
        data = [[ctree.compare(otree, method='RF'), ctree.compare(otree, method='MRCA')] for otree in trimmed_forest]
        rf_data[extra_name + kwargs.name] = [l[0] for l in data]
        mrca_data[extra_name + kwargs.name] = [l[1] for l in data]

    # TODO also check single trees found by cli runs of gctree

    with open(base_inference_name + '.inference.1.p', 'rb') as fh:
        base_ctree = pickle.load(fh)
    rf_data['cli_base'] = [ctree.compare(base_ctree, method='RF')]
    mrca_data['cli_base'] = [ctree.compare(base_ctree, method='MRCA')]
    with open(mut_inference_name + '.inference.1.p', 'rb') as fh:
        fancy_ctree = pickle.load(fh)
    rf_data['cli_all_criteria'] = [ctree.compare(fancy_ctree, method='RF')]
    mrca_data['cli_all_criteria'] = [ctree.compare(fancy_ctree, method='MRCA')]
    with open(output_dir / 'mrca_data.p', 'wb') as fh:
        fh.write(pickle.dumps(mrca_data))

    with open(output_dir / 'rf_data.p', 'wb') as fh:
        fh.write(pickle.dumps(rf_data))

    print("==========================")
    print(rf_data)
    print(fancy_ctree.tree)
    print(ctree.tree)
        

if __name__ == "__main__":
    main()
