
import argparse
import ete3 as ete 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse():
    parser = argparse.ArgumentParser("extract nucleotide FASTA files from Panaroo gene_data.csv output")
    parser.add_argument("gene_data", help="Panaroo output")
    parser.add_argument("fasta_outdir", help="nucleotide FASTA")
    parser.add_argument("--ids", help="restrict output to ids")
    args = parser.parse_args()
    extract_sequences(**vars(args))

def extract_sequences(gene_data, fasta_outdir, ids):
    if ids:
        with open(ids, 'r') as f:
            assembly_ids = []
            for l in f:
                assembly_id = l.strip()
                assembly_ids.append(assembly_id)
        assembly_ids = set(assembly_ids)
    with open(gene_data, 'r') as f:
        f.readline()
        current_out = None
        current_gff = None
        for l in f:
            try:
                gff_file, scaffold_name, clustering_id, annotation_id, prot_sequence, dna_sequence, gene_name, description = l.strip().split(",")
            except ValueError:
                print(l)
            if not ids or gff_file in assembly_ids:
                if current_out is None:
                    current_out = open("{}/{}.fna".format(fasta_outdir, gff_file), 'a')
                elif current_gff != gff_file:
                    current_out = open("{}/{}.fna".format(fasta_outdir, gff_file), 'a')
                current_gff = gff_file 
                seq = Seq(dna_sequence)
                record = SeqRecord(seq, annotation_id, description = description)
                SeqIO.write(record, current_out, "fasta")

if __name__ == "__main__":
    parse()
