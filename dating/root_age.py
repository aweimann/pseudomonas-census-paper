import argparse
import os
import dendropy as dp
import Bio.AlignIO
import pandas as pd

def parse():
    parser = argparse.ArgumentParser("root age")
    parser.add_argument("real_tree", help="treeannotator generated tree based on real dates")
    parser.add_argument("dates",  help = "dates")
    parser.add_argument("alignment",  help = "alignment")
    parser.add_argument("outdir", help="outdir for results")
    parser.add_argument("st", help="clone name")
    args = parser.parse_args()
    process(**vars(args))

def process(real_tree, dates, alignment, st, outdir):
    #read in alignment and get alignment length
    aln = Bio.AlignIO.read(alignment, format = "fasta")
    aln_length = aln.get_alignment_length()
    samples = [i.name for i in aln._records]
    dates = pd.read_csv(dates, sep = "\t", index_col = 0, parse_dates = ['collection_date'])
    dates.collection_date = dates.collection_date.apply(lambda x: x.year + x.dayofyear/365)
    dates = dates.loc[samples, ]
    dates.sort_values(by = "collection_date", inplace =True, ascending = False)
    youngest_age = dates.iloc[0].loc["collection_date"]
    #read in trees and extract root age
    root_ages = []
    root_ages_hpd = []
    real_root_age = None
    real_root_age_hpd = None
    with open(real_tree) as f:
        tree = dp.Tree.get(file = f, schema = "nexus")
        new, old  = tree.seed_node.annotations["height_95%_HPD"]._get_value()
        real_root_age_hpd = (youngest_age - float(new), youngest_age - float(old))
        real_root_age = youngest_age - float(tree.seed_node.annotations["height_median"]._get_value())
    real_root_age_out = "%s/real_root_age.txt" % outdir
    #read in dates and get age of youngest tip
    if os.path.exists(real_root_age_out):
        f = open(real_root_age_out, 'a')
        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (st, real_root_age, real_root_age_hpd[0], real_root_age_hpd[1], aln_length, youngest_age))
    else:
        f = open(real_root_age_out, 'w')
        f.write("st\troot_age\troot_age_2.5perc\troot_age_97.5_perc\talignment_length\tyoungest_age\n")
        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (st, real_root_age, real_root_age_hpd[0], real_root_age_hpd[1], aln_length, youngest_age))
    root_ages_out = "%s/root_ages.txt" % outdir
    f.close()


if __name__ == "__main__":
    parse()
