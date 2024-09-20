import click
import json
from collections import Counter
import gctree
import pickle
import historydag as hdag
from pathlib import Path
import ete3

# mutability file
# substitution file
# parsimony forest
DEBUG=False


def compute_reversions_of_history(history):
    reversions = 0
    edges = list(history.get_edges(skip_ua_node=True))
    root = edges[0][0]
    root.left_bases = [frozenset()] * len(root.label.sequence)
    for parent, child in edges:
        child_left_bases = parent.left_bases.copy()
        for idx, (pnuc, cnuc) in enumerate(zip(parent.label.sequence, child.label.sequence)):
            if pnuc != cnuc:
                if cnuc in parent.left_bases[idx]:
                    reversions += 1
                child_left_bases[idx] |= frozenset({pnuc})
        child.left_bases = child_left_bases
    return reversions


@click.command()
@click.argument('output_file')
@click.argument('parsimony_forest')
@click.argument('fivemer_mutabilities')
@click.argument('fivemer_substitution')
@click.argument('input_sequences_path')
@click.argument('meta_file')
@click.argument('true_treespath')
@click.argument('dnapars_outfile')
@click.argument('abundances')
@click.argument('root_name')
@click.option('-a', '--all_dagtrees_data', default=None, type=str)
def main(output_file, parsimony_forest, fivemer_mutabilities, fivemer_substitution, input_sequences_path, meta_file, true_treespath, dnapars_outfile, abundances, root_name, all_dagtrees_data):
    parsimony_forest = Path(parsimony_forest)
    output_file = Path(output_file)
    with open(parsimony_forest, 'rb') as fh:
        forest = pickle.load(fh)
    dag = forest._forest
    dnapars_trees = gctree.phylip_parse.parse_outfile(dnapars_outfile, abundances, root_name)
    naive_seq = forest._validation_stats["root_seq"]

    with open(meta_file, 'r') as fh:
        chain_split = json.loads(fh.read())["l_offset"]

    ll_dagfuncs = gctree.branching_processes._ll_genotype_dagfuncs(*forest.parameters).weight_funcs
    mut_funcs = gctree.branching_processes._mutability_dagfuncs(
        splits=[int(chain_split)],
        mutability_file=fivemer_mutabilities,
        substitution_file=fivemer_substitution,
    ).weight_funcs
    poisson_funcs = gctree.mutation_model._context_poisson_likelihood_dagfuncs(
        splits=[int(chain_split)],
        mutability_file=fivemer_mutabilities,
        substitution_file=fivemer_substitution,
    ).weight_funcs
    allele_funcs = hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": lambda n1, n2: 0 if n1.is_ua_node() else int(n1.label.sequence != n2.label.sequence),
            "accum_func": sum,
        },
        name="NumAlleles",
    )
    reversion_funcs = gctree.branching_processes._naive_reversion_dagfuncs(naive_seq).weight_funcs
    placeholder_funcs = hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": lambda n1, n2: 0,
            "accum_func": sum,
        },
        name="Whole DAG",
    )
    
    kwargls = [
        ll_dagfuncs,
        poisson_funcs,
        combined_funcs := ll_dagfuncs + poisson_funcs,
        rev_combined_funs := reversion_funcs + ll_dagfuncs + poisson_funcs,
        no_bp_funcs := reversion_funcs + poisson_funcs,
        placeholder_funcs,
    ]
    combined_funcs.name = "Likelihood_then_Context"
    rev_combined_funs.name = "Reversions_then_Default"
    no_bp_funcs.name = "Reversions_then_Context"



    try:
        modeltree, matched_simu_path = get_true_tree(input_sequences_path, true_treespath)
    except Exception as e:
        print("Error finding true tree", input_sequences_path, e)
        if DEBUG:
            raise e
        return

    # modeltree.summary()
    # dag.summary()
    # (modeltree | dag).summary()
    ts1 = {n.label for n in modeltree.get_leaves()}
    ts2 = {n.label for n in dag.get_leaves()}
    if  ts1 != ts2:
        # This is sketchy because I don't understand why it happens. May affect
        # parsimony score but shouldn't affect RF distance (which is only thing
        # true tree is used for!)

        t1_err = list(ts1 - ts2)
        t2_err = list(ts2 - ts1)
        if len(t1_err) == len(t2_err) and len(t1_err) == 1:
            old_label = t1_err[0]
            new_label = t2_err[0]
            # This guarantees the nodes are actually matched
            assert new_label.abundance == 0

            def l_func(in_node):
                in_label = in_node.label
                if in_label == old_label:
                    return new_label
                else:
                    return in_label

            modeltree = modeltree.relabel(l_func)
            ts1 = {n.label for n in modeltree.get_leaves()}
            if ts1 != ts2:
                print("Error encountered (maybe convergent mutation)", parsimony_forest)
                print(ts1 - ts2)
                print("=======")
                print(ts2 - ts1)
                if DEBUG:

                    raise RuntimeError
                return

        else:
            print("Error encountered (maybe convergent mutation)", parsimony_forest)
            print(ts1 - ts2)
            print("=======")
            print(ts2 - ts1)
            if DEBUG:

                raise RuntimeError
            return

    
    node_count_funcs = hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": lambda n1, n2: 0 if n1.is_ua_node() else 1,
            "accum_func": sum,
        },
        name="NumNodes",
    )
    true_tree_comparison_funcs = (hdag.utils.make_rfdistance_countfuncs(modeltree, rooted=True) + 
                                  node_count_funcs)
    true_tree_comparison_funcs.name = "RootedRF_and_nodecount"

    tree_summary_funcs = true_tree_comparison_funcs + ll_dagfuncs + mut_funcs + allele_funcs + poisson_funcs
    tree_summary_col_names = tree_summary_funcs.names

    # Throw out non-unique trees...
    dnapars_histories = []
    dnapars_unique_set = set()
    for dptree in dnapars_trees:
        dpforest = gctree.branching_processes.CollapsedForest([gctree.phylip_parse.disambiguate(dptree)])
        node_set = frozenset(dpforest._forest.preorder(skip_ua_node=True))
        if node_set not in dnapars_unique_set:
            dnapars_unique_set.add(node_set)
            dnapars_histories.append(dpforest._forest)

    dnapars_data = [dphistory.optimal_weight_annotate(**tree_summary_funcs) for dphistory in dnapars_histories]

    rf_data = {}
    dag_weight_values = {}
    for kwargs in kwargls:
        trimmed_dag = dag.copy()
        weight_val = trimmed_dag.trim_optimal_weight(**kwargs, optimal_func=min)
        data = trimmed_dag.weight_count(**true_tree_comparison_funcs)
        rf_data[kwargs.name] = data
        dag_weight_values["best_" + kwargs.name] = weight_val

    # Get reversions for default ranking
    rev_dag = dag.copy()
    rev_dag.trim_optimal_weight(**combined_funcs)
    rev_dag = rev_dag[0]

    with open(output_file, 'wb') as fh:
        fh.write(pickle.dumps(
            {"WholeDAG": rf_data,
             "WholeDAGTrimVals": dag_weight_values,
             "NumLeaves": modeltree.num_leaves(),
             "TrueTreeNumNodes": modeltree.optimal_weight_annotate(**node_count_funcs),
             "dnaparsTrees": [tree_summary_col_names] + dnapars_data,
             "SimulationPath": matched_simu_path,
             "InferencePath": parsimony_forest,
             "LikelihoodThenContextReversions": compute_reversions_of_history(rev_dag),
             "SimuReversions": compute_reversions_of_history(modeltree)}
        ))


    n_trees = dag.count_histories()
    if all_dagtrees_data is not None:
        dag_trees_data = [tree_summary_col_names] + [history.optimal_weight_annotate(**tree_summary_funcs) for history in dag]
        with open(all_dagtrees_data, 'wb') as fh:
            fh.write(pickle.dumps(
                {
                    "dnaparsTrees": [tree_summary_col_names] + dnapars_data,
                    "dagtrees": dag_trees_data,
                    "NumLeaves": modeltree.num_leaves(),
                    "TrueTreeNumNodes": modeltree.optimal_weight_annotate(**node_count_funcs),
                }
            ))


    # print("==========================")
    if DEBUG:
        print(rf_data)
        name_dict = {leaf: str(idx) for idx, leaf in enumerate(dag.get_leaves())}
        hdag.dag.ascii_compare_histories(dag[0], modeltree, lambda n: '' if n not in name_dict else name_dict[n], sort_method="leaf-name")
        # ladderize, leaf-name, child-name
        # compact=True
        

def read_true_newick(newicks_path, match_path):
    with open(newicks_path, 'r') as fh:
        for line in fh:
            path, newick = line.split(' ')
            if Path(path).resolve().samefile(match_path.resolve()):
                return newick
    raise RuntimeError(f"No newick matching path {match_path} found in file {newicks_path}.")


def convert_inference_seqname(key):
    # examples: 
    # leaf-abcdefhijlpqrstuvxyz_contig_igh+igk   --> leaf-abcdefhijlpqrstuvxyz
    # leaf-abcdefhiklmnoprstuwx-2_contig_igh+igk --> leaf-abcdefhiklmnoprstuwx
    parts = key.split('_')[0].split('-')
    return parts[0] + '-' + parts[1]


def get_true_tree(input_sequences_path, true_treespath):
    input_sequences_path = Path(input_sequences_path)
    input_fasta = hdag.utils.load_fasta(input_sequences_path)
    simu_path = input_sequences_path.parent / "../../../simu/selection/simu/"
    possible_paths = list(simu_path.glob("event-*"))

    def evaluate_paths(possible_paths):
        for path in possible_paths:
            candidate_fasta = hdag.utils.load_fasta(path / "simu.fasta")
            if all(convert_inference_seqname(key) in candidate_fasta for key in input_fasta if key != "XnaiveX"):
                return path
        raise RuntimeError("Could not match simulation with inference")

        
    matched_path = evaluate_paths(possible_paths)


    simu_tree = ete3.Tree(newick=read_true_newick(true_treespath, matched_path), format=1)


    # Double check that converted sequences are the same for leaves (also need
    # to convert leaf name keys)

    # get conversion function because sometimes padding Ns are mutated
    firstn_idx = simu_tree.nuc_seq.find('N')
    lastn_idx = simu_tree.nuc_seq.rfind('N')
    if firstn_idx == -1:
        def seq_convert(seq):
            return seq
    else:
        def seq_convert(seq):
            return seq[:firstn_idx] + seq[lastn_idx + 1:]

    # put sequences from simu_alignment on sim_tree
    for node in simu_tree.traverse():
        node.add_feature("sequence", seq_convert(node.nuc_seq))
    s1 = {n.name: n.sequence for n in simu_tree.iter_leaves()} 
    s2 = {convert_inference_seqname(key): val for key, val in input_fasta.items() if key != 'XnaiveX'}
    if s1 != s2:
        print("keys match", set(s1.keys()) == set(s2.keys()))
        print("seqs match", set(s1.values()) == set(s2.values()))
        print("num seqs", len(s1), len(s2))
        print("num unique seqs", len(set(s1.values())), len(set(s2.values())))
        print("seq lengths:", len(next(iter(s1.values()))), len(next(iter(s2.values()))))
        set1 = set(s1.keys())
        set2 = set(s2.keys())
        print(set1 - set2, set2 - set1)
        raise RuntimeError("Couldn't find tree with matching leaf sequences")

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
    return (cforest._forest, matched_path)

if __name__ == "__main__":
    main()
