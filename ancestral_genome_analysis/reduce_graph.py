import argparse
import networkx as nx

def parse():

    parser = argparse.ArgumentParser("reduce Panaroo pan-genome graph")
    parser.add_argument("graph", help="Panaroo graph")
    parser.add_argument("isolates", help="isolates to keep")
    parser.add_argument("out_graph", help="output reduced graph")
    args = parser.parse_args()
    reduce(**vars(args))

def reduce(graph, isolates, out_graph):
    g = nx.read_gml(graph)
    g_new = nx.Graph()
    isols = set()
    with open(isolates) as f:
        #skip header
        f.readline()
        for l in f:
            isols.add(l.strip())
    #isolate id to name dict
    id2isol = {}
    for isol, i in zip(g.graph["isolateNames"], range(len(g.graph["isolateNames"]))):
            id2isol[str(i)] = isol
    isol2id = dict([(j, i) for i, j in id2isol.items()])
    for isol in isols:
        if isol not in isol2id:
            print(isol, " missing")
    isols = list(filter(lambda isol: isol in isol2id.keys(), isols))
    g_new.graph["isolateNames"] = [isol for isol in g.graph["isolateNames"] if isol in isols]
    for n in g.nodes:
        if not "members" in g.nodes[n]:
            g_new.add_node(n, **node_attrs)
            continue

        if type(g.nodes[n]["members"]) == type([]):
            members = g.nodes[n]["members"]
            seqIDs = g.nodes[n]["seqIDs"]
            lengths = g.nodes[n]["lengths"]

        else:
            members = [g.nodes[n]["members"]]
            seqIDs = [g.nodes[n]["seqIDs"]]
            lengths = [g.nodes[n]["lengths"]]
        members_filt = [j for j in filter(lambda i: id2isol[str(i)] in isols, members)]
        seqIDs_filt = [j for j in filter(lambda i: id2isol[i.split("_")[0]] in isols, seqIDs)]
        geneIDs = [ j for j in filter(lambda i: id2isol[i.split("_")[0]] in isols, g.nodes[n]["geneIDs"].split(";"))]
        genomeIDs = [j for j in filter(lambda i: id2isol[str(i)] in isols, g.nodes[n]["genomeIDs"].split(";"))]
        lengths_indices =  [j for j in filter((lambda i: id2isol[str(i[0])] in isols), zip(members, range(len(members))))]
        if len(list(members_filt)) != 0:
            node_attrs = g.nodes[n]
            node_attrs["members"] = list(members_filt)
            node_attrs["seqIDs"] = list(seqIDs_filt)
            node_attrs["genomeIDs"] = ";".join(list(genomeIDs))
            node_attrs["geneIDs"] = ";".join(list(geneIDs))
            node_attrs["lengths"] = [] 
            node_attrs["size"] = len(node_attrs["members"])
            g_new.add_node(n, **node_attrs)
    for e in g.edges:
        if not "members" in g.edges[e]:
            g_new.add_edge(e[0], e[1])
            continue
        if e[0] in g_new.nodes and e[1] in g_new.nodes:
            if type(g.edges[e]["members"]) == type([]):
                members = g.edges[e]["members"]
            else:
                members = [g.edges[e]["members"]]
            members_filt = [j for j in filter(lambda i: id2isol[str(i)] in isols, members)]
            genomeIDs = [j for j in filter(lambda i: id2isol[str(i)] in isols, g.edges[e]["genomeIDs"].split(";"))]
            if len(list(members_filt)) != 0:
                edge_attrs = g.edges[e]
                edge_attrs["members"] = list(members_filt)
                edge_attrs["genomeIDs"] = list(genomeIDs)
                edge_attrs["weight"] = len(edge_attrs["members"])
                g_new.add_edge(e[0], e[1], **edge_attrs)
        
    nx.readwrite.gml.write_gml(g_new,out_graph, stringizer = str )
    
    
if __name__ == "__main__":
    parse()
