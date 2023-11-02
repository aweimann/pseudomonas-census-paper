import pandas as pd
import argparse
import sys
sys.setrecursionlimit(10000)
import ete3
sys.path.append('/lustre/scratch118/infgen/team216/aw27/pa_clinical/results/parsimony/')
import parsimony



def parse():
    parser = argparse.ArgumentParser("Panaroo stats across all clones")
    parser.add_argument("gpas", help="gene presence/absence Rtab")
    parser.add_argument("isol_info", help="assembly info")
    parser.add_argument("st", help="clone info")
    parser.add_argument("tree", help="clone type phylogenetic tree")
    parser.add_argument("outdir", help="summary output table")
    args = parser.parse_args()
    get_closest_representative(**vars(args))



def prune_tree(tree, isol_info, st, outdir):
    tree = ete3.Tree(tree, format = 1)
    isols = []
    isols = pd.read_csv(isol_info, sep = "\t", header = 0)
    isols = isols.loc[isols.loc[:, "majority_ST"] == st, ]
    isols_in_tree = [] 
    isols = set(isols.loc[:, "sample_name"])
    for leaf in tree.iter_leaves():
        if leaf.name in isols:
            isols_in_tree.append(leaf.name)
    tree.prune(isols_in_tree)
    tree.write(outfile = outdir + "/{}_pruned_tree.nwk".format(st), format = 1)
    return isols

def subsample_panaroo_rtab(rtab, isol_info, outdir, st, isols_in_tree):
    t = pd.read_csv(rtab, sep = "\t", index_col = 0)
    t.columns = [i.replace(".spades", "") for i in t.columns]
    names = pd.read_csv(isol_info, sep = "\t", header = 0)
    isols = names.loc[:, "gff"].apply(lambda x: x.split("/")[-1].split(".gff")[0])
    names = names.assign(panaroo_id = isols)
    names = names.loc[names.loc[:, "majority_ST"] == st, :]
    names_filter  = [i in isols_in_tree for i in names.loc[:, "sample_name"]]
    names = names.loc[names_filter, ]
    t = t.loc[:, names.panaroo_id]
    t.columns = names.loc[:, "sample_name"]
    t.T.to_csv(outdir + "/{}_gene_presence_absence_transposed.Rtab".format(st), sep = "\t")

def get_closest_representative(gpas, isol_info, tree, st, outdir):
    isols_in_tree = prune_tree(tree, isol_info, st, outdir)
    subsample_panaroo_rtab(gpas, isol_info, outdir, st, isols_in_tree)
    parsimony.parsimony(outdir + "/{}_pruned_tree.nwk".format(st), outdir + "/{}_gene_presence_absence_transposed.Rtab".format(st), outdir + "/parsimony_{}/".format(st), tree_is_named = True, is_transposed = False)
    anc = pd.read_csv(outdir + "/parsimony_{}/ancestral_states.txt".format(st), sep = "\t", header = 0, index_col = 0)
    anc_sim = anc.apply(lambda i: ((anc.iloc[0,:] == 1) & (i == 1)).sum(), axis = 1).sort_values()
    anc_sim_terminal = anc_sim.loc[[not "internal" in str(i) and not pd.isnull(i) for i in anc_sim.index]]
    anc_sim_max = anc_sim_terminal.loc[anc_sim_terminal == anc_sim_terminal.max()]
    anc_sim_max = pd.DataFrame(anc_sim_max.iloc[0:1, ])
    anc_sim_max = anc_sim_max.assign(anc = anc_sim.max(), majority_ST = st)
    anc_sim_max.to_csv(outdir + "/closest_isolate.txt".format(st), header = False, mode = "a", sep = "\t")

if __name__ == "__main__":
    parse()
