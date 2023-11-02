import argparse
import pandas as pd
from Bio import SeqIO, Seq
import os

def parse():
    parser = argparse.ArgumentParser("subset alignmenment based on an input table of samples and sequence types")
    parser.add_argument("out_dir", help="output directory")
    parser.add_argument("sample_annotation", help="sequence type annotation table")
    parser.add_argument("group", help="use the column to group alignment")
    parser.add_argument("--reference", help="reference fasta")
    parser.add_argument("--alns", nargs = "*", help="pseudoalignments")
    parser.add_argument("--single_output", help="name of single output fasta file")
    args = parser.parse_args()
    subsample(**vars(args))

def subsample(out_dir, reference, sample_annotation, alns, single_output, group):
    annot = pd.read_csv(sample_annotation, sep = "\t", index_col = 0, header = 0, engine = 'python')
    annot = annot.astype({group: "str"})
    st_seen = {}
    if single_output:
        if reference:
            with open(reference, 'r') as ref:
                record_iter_ref = SeqIO.parse(ref, "fasta")
                for entry_ref in record_iter_ref:
                    entry_ref.id = "PAO1"
                with open(single_output, 'w') as f:
                    SeqIO.write(entry_ref, f, "fasta")
    for aln in alns:
        with open(aln, 'r') as aln:
            record_iter = SeqIO.parse(aln, "fasta")
            for entry in record_iter:
                entry_cleaned = entry.name.replace("_AE004091","")
                if entry_cleaned in annot.index:
                    ST = annot.loc[entry_cleaned, group]
                    sample_name = annot.loc[entry_cleaned, "sample_name"]
                    print(ST, sample_name)
                    out = entry
                    entry.id = sample_name
                    if single_output:
                        with open(single_output, 'a') as f:
                            SeqIO.write(entry, f, "fasta")
                    else:
                        out_fasta ="%s/%s.fasta"  %(out_dir, ST)
                        if not ST in st_seen :
                            with open(out_fasta, 'w') as f:
                                st_seen[ST] = True
                                if reference:
                                    with open(reference, 'r') as ref:
                                        record_iter_ref = SeqIO.parse(ref, "fasta")
                                        for entry_ref in record_iter_ref:
                                            entry_ref.id = "PAO1"
                                        SeqIO.write(entry_ref, f, "fasta")
                        with open(out_fasta, 'a') as f:
                            SeqIO.write(entry, f, "fasta")

if __name__ == "__main__":
    parse()
