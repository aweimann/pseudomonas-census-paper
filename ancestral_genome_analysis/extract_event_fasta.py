import pandas as pd
import networkx as nx
import argparse
import sys

def parse():
    parser = argparse.ArgumentParser("Panaroo stats across all clones")
    parser.add_argument("events", help="summary with gene gain/loss events")
    parser.add_argument("panaroo_graph", help="Panaroo final graph")
    parser.add_argument("out_fasta", help="summary output table")
    parser.add_argument("--split_events", action = "store_true", help="split events into individual genes")
    args = parser.parse_args()
    extract_fasta(**vars(args))


def extract_fasta(events, panaroo_graph, out_fasta, split_events):
    events = pd.read_csv(events, sep = "\t", header= 0)
    gene2events = {}
    for index, row in events.iterrows():
        if row.loc["Number.of.genes"] > 1:
            genes = row.loc["Genes"].split(',')
            event = row.loc["Event"]
            for gene in genes:
                if gene not in gene2events:
                    gene2events[gene] = [event]
                else:
                    gene2events[gene].append(event)
    
    event2genes = {}
    for gene, events in gene2events.items():
        for event in events:
            if event not in event2genes:
                event2genes[event] = [gene]
            else:
                event2genes[event].append(gene)

    gene2dna = {}
            
    g = nx.read_gml(panaroo_graph)
    for n in g.nodes:
        gene_name = g.nodes[n]["name"]
        if gene_name in gene2events:
            if split_events:
                dna = g.nodes[n]["protein"].split(';')[0]
            else:
                dna = g.nodes[n]["dna"].split(';')[0]
            gene2dna[gene_name] = dna
    with open(out_fasta, 'w') as out_fasta:
        if split_events:
            for event, genes in event2genes.items(): 
                    for gene in genes:
                        out_fasta.write(">{}_{}\n".format(event, gene))
                        out_fasta.write("{}\n".format(gene2dna[gene]))

        else:
                for event, genes in event2genes.items(): 
                    if len(genes) > 5:
                        out_fasta.write(">{}\n".format(event))
                        for gene in genes:
                            out_fasta.write("{}".format(gene2dna[gene]))
                        out_fasta.write("\n")


if __name__ == "__main__":
    parse()
