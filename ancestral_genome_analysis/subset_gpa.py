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
    parser.add_argument("tree", help="clone type phylogenetic tree")
    parser.add_argument("out_table", help="summary output table")
    args = parser.parse_args()
    subsample_panaroo_rtab(**vars(args))


def subsample_panaroo_rtab(gpas, isol_info, tree, out_table):
    tree = ete3.Tree(tree, format = 1)
    tree_isols = [n.name for n in tree.iter_leaves()]
    names = pd.read_csv(isol_info, sep = "\t", header = 0)
    isols = names.loc[:, "gff"].apply(lambda x: x.split("/")[-1].split(".gff")[0])
    names = names.assign(panaroo_id = isols)

    names_filter  = [i in tree_isols for i in names.loc[:, "sample_name"]]
    names = names.loc[names_filter, ]
    t = pd.read_csv(gpas, sep = "\t", index_col = 0)
    t.columns = [i.replace(".spades", "") for i in t.columns]
    t = t.loc[:, names.panaroo_id]
    t.columns = names.loc[:, "sample_name"]
    t.T.to_csv(out_table, sep = "\t")

if __name__ == "__main__":
    parse()
