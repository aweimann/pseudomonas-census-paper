import argparse
import dendropy as dp
import random

def parse():
    parser = argparse.ArgumentParser("process table of pathoadaptive mutations and generate samples overview")
    parser.add_argument("patho_table_in",  help="table with pathoadaptive mutations per node and sequence cluster")
    parser.add_argument("tree_dir",  help="directory with a phylogenetic per sequence cluster")
    parser.add_argument("patho_table_out", help="output table sample by gene with all pathoadaptive mutations found in each samples")
    parser.add_argument("sample_trajectories_out", help="individual trajectories per sample")
    parser.add_argument("node2depth_out", help="output table sample by gene with all pathoadaptive mutations found in each samples")
    parser.add_argument("gene2position_out", help="output table sample by gene with all pathoadaptive mutations found in each samples")
    args = parser.parse_args()
    process(**vars(args))

def node_split(node_st):
    if node_st is None:
        return None, None
    node_sep = node_st.split('_')
    node = "_".join(node_sep[:-1])
    st = node_sep[-1]
    return node, st

def process(patho_table_in, tree_dir, patho_table_out, sample_trajectories_out, node2depth_out, gene2position_out):
    #keep track of the position in pathoadpative trajectory for each node
    gene2position = {}
    node2position = {}

    node_st2patho = {}
    sts = set()
    genes = set()
    with open(patho_table_in) as f:
        f.readline()
        #iterate through pathoadaptive mutations and create dictonary of node, st to gene
        for line in f:
            node, PAO1, gene_name, st = line.strip().split("\t") 
            sts.add(st)
            genes.add("%s_%s" % (PAO1, gene_name))
            node_st = "%s_%s" % (node, st)
            if node_st not in node_st2patho:
                node_st2patho[node_st] = ["%s_%s" % (PAO1, gene_name)]
            else:
                node_st2patho[node_st].append("%s_%s" % (PAO1, gene_name))
    sample2patho = {}
    patho_table_out =  open(patho_table_out, 'w')
    node2depth = {}
    sample_trajectories_out = open(sample_trajectories_out, 'w')
    #iterate over all clades/STs
    for st in sts:
        #assign depth to each node
        tree = dp.Tree.get(path = "%s/%s_rescaled.nwk" % (tree_dir, st), schema = "newick")
        for node in  tree.preorder_node_iter():
            if node.parent_node is None:
                node.label = "root"
                label_st = node.label + "_" + st
                sample2patho[label_st] = set()
                node2position[label_st] = 0
                node2depth[label_st] = (0, None)
                continue
            if node.taxon is None:
                node_name = node.label.replace(" ", "_")
            else:
                node_name = node.taxon.label.replace(" ", "_")

            node_st = "%s_%s" % (node_name, st)
            if node_st in node_st2patho:
                node_set = set(node_st2patho[node_st])
            else:
                node_set = set()
            out = []
            for gene in node_set:
                out.append(gene)
            parent_node = node.parent_node
            label = parent_node.label.replace(" ", "_")
            label_st = "%s_%s" % (label, st)
            node2depth[node_st] = (node2depth[label_st][0] + 1, label_st)
            parent_node_set = set(sample2patho[label_st])
            sample2patho[node_st] =  node_set.union(parent_node_set)
            #update node position in trajectory randomly for mutations happening indistinguishable at the same time/position
            # keep node, st information for later filtering
            for gene in node_set:
                if gene not in gene2position:
                    gene2position[gene] = [(node2position[label_st] + random.randint(1, len(node_set)), node_name, st)]
                else:
                    gene2position[gene].append((node2position[label_st] + random.randint(1, len(node_set)), node_name, st))
            node2position[node_st] = node2position[label_st] + len(node_set)
        for tip in tree.leaf_node_iter():
            out = []
            tip_name = tip.taxon.label.replace(" ", "_")
            tip_st = "%s_%s" % (tip_name, st)
            if tip_st in node_st2patho:
                node_set = set(node_st2patho[tip_st])
            else:
                node_set = set()
            sample_trajectories_out.write(tip.taxon.label + "\t")
            sample_trajectories_out.write(",".join(node_set) + "-")
            current_node = tip.parent_node
            while not "root" in current_node.label:
                label = current_node.label.replace(" ", "_")
                label_st = "%s_%s" % (label, st)
                out = []
                if label_st in node_st2patho:
                    node_set = set(node_st2patho[label_st])
                else:
                    node_set = {}
                for gene in node_set:
                    out.append(gene)
                sample_trajectories_out.write(",".join(out) + "-")
                current_node = current_node.parent_node
            sample_trajectories_out.write("\n")
    sample_trajectories_out.close()
    patho_table_out.write("sample_name\tmajority_ST\t")
    patho_table_out.write("\t".join(list(genes)) + "\n")
    for sample_st in sample2patho:
        sample, st = node_split(sample_st) 
        patho_table_out.write(sample + "\t" + st +  "\t")
        out = [] 
        for gene in genes:
            if gene in sample2patho[sample_st]:
                out.append("1")
            else:
                out.append("0")
        patho_table_out.write("\t".join(out) + "\n")
    patho_table_out.close()
    with open(node2depth_out, 'w') as f:
        for node_st in node2depth:
            node,st = node_split(node_st)
            depth, ancestor_st = node2depth[node_st]  
            ancestor, st = node_split(ancestor_st)
            if ancestor is None:
                continue
            f.write("{}\t{}\t{}\t{}\n".format(ancestor, node, st, depth))
    with open(gene2position_out, 'w') as f:
        f.write("gene\tposition\node\st\n")
        for gene, positions in gene2position.items():
            for position in positions:
                f.write("{}\t{}\t{}\t{}\n".format(gene, position[0], position[1], position[2]))


if __name__ == "__main__":
    parse()
