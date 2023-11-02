import pandas as pd
import argparse
import sys

def parse():
    parser = argparse.ArgumentParser("parse emapper output and assign to sts")
    parser.add_argument("emapper", help="eggNOG mapper output")
    parser.add_argument("event_summary", help = "event_summary")
    parser.add_argument("clone2function", help="summary output table")
    parser.add_argument("annotation_type", choices = ["KEGG_Pathway", "COG_category"], help = "event_summary")
    args = parser.parse_args()
    map_clones(**vars(args))

def map_clones(emapper, event_summary, clone2function, annotation_type):
    event2function = {}
    event2clone = {}
    event2event_type = {}
    with open(emapper) as f:
        #skip header
        for i in range(4):
            f.readline()
        header = f.readline().strip().split('\t')
        field_names = dict([(i, j) for i, j in zip(header, range(len(header)))])
        print(field_names)
        for l in f:
            if l.startswith("#"):
                continue
            fields  = l.strip().split('\t')
            COGs = fields[field_names[annotation_type]]
            if annotation_type == "KEGG_Pathway":
                COGs = COGs.split(",")
            event = fields[field_names["#query"]].split("_")[0] 
            gene = "_".join(fields[field_names["#query"]].split("_")[1:])
            if COGs[0] == "-" or len(COGs) == 0:
                continue
            else:
                for COG in COGs:
                    if event in event2function:
                        event2function[event].append((COG, gene))
                    else:
                        event2function[event] = [(COG, gene)]
    with open(event_summary) as f:
        header = f.readline().strip().split('\t')
        field_names = dict([(i, j) for i, j in zip(header, range(len(header)))])
        for l in f:
            fields  = l.strip().split('\t')
            event = fields[field_names["Event"]]
            upstream, downstream = fields[field_names["Branch"]].split(":")
            gain_loss = fields[field_names["Gain.or.loss"]]
            event2clone[event] = downstream
            event2event_type[event] =gain_loss 
    functions = set()
    with open(clone2function, 'w') as out:
        out.write("clone\tevent\tevent_type\tCOG_category\tgene_family\n")
        for event, functions in event2function.items():
            clone = event2clone[event]
            event_type = event2event_type[event]
            for function, gene in functions:
                out.write("{}\t{}\t{}\t{}\t{}\n".format(clone, event, event_type, function, gene))








if __name__ == "__main__":
    parse()
