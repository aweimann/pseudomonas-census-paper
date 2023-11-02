
import argparse
import dendropy as dp

def parse():
    parser = argparse.ArgumentParser("stratify branches by transmission type")
    parser.add_argument("sample2patient", help="sample to patient")
    parser.add_argument("sample2infection", help="sample to infection type")
    parser.add_argument("branches_out", help="classification of branches")
    parser.add_argument("trees", nargs = "*",  help="one or more rooted phylogenetic trees")
    args = parser.parse_args()
    extract(**vars(args))

def extract(sample2patient, sample2infection, trees, branches_out):
    sample2patient_dict = {}
    patient2samples = []
    sample2infection_dict = {}
    patient2st = {}
    with open(sample2patient) as f:
        for l in f:
            sample, st, patient = l.strip().split('\t')
            patient2st[patient] = st
            sample2patient_dict[sample]  = patient
    with open(sample2infection) as f:
        for l in f:
            sample, infection = l.strip().split('\t')
            sample2infection_dict[sample]  =infection 
    branches = []
    #output file
    out_type = open(branches_out, 'w')
    out_type.write("node\tlevel7000\tinfection_type\tis_within_patient\tis_cf_transmissible\tis_non_cf_transmissible\n")
    #get all patients in this tree
    for tree in trees:
        node2patients_dict = {}
        dp_tree = dp.Tree.get(path = tree, schema = 'newick', preserve_underscores = True)
        #iterate over nodes in order of distance from root
        levelorder = [i for i in dp_tree.levelorder_node_iter()]
        levelorder.reverse()
        #assign all patients to node that are downstream
        for node in levelorder:
            if not node.taxon is None:
                if node.taxon.label == "PAO1": 
                    continue
                name = node.taxon.label
                patient = sample2patient_dict[name]
                node2patients_dict[name] = [patient]
            else:
                #get all patients from nodes downstream
                name = node._label
                patients = set()
                for leaf in node.leaf_iter():
                    if leaf.taxon.label == "PAO1":
                        continue
                    patients.add(sample2patient_dict[leaf.taxon.label])
                node2patients_dict[name] = list(patients)
        for node, patients in node2patients_dict.items():
            infection_types = set()
            for patient in patients:
                infection_types.add(sample2infection_dict[patient])
            is_within_patient = len(patients) == 1
            is_cf_transmissible = "cystic fibrosis" in infection_types and not is_within_patient
            #handle NAs appropriately
            is_non_cf_transmissible = False
            infection_type = None
            if len(infection_types) > 1:
                infection_type = "multiple"
            else:
                infection_type = infection_types.pop()
                infection_types.add(infection_type)
            if "cystic fibrosis" in infection_types:
                if "NA" in infection_types:
                    if len(infection_types) > 2:
                        is_non_cf_transmissible = True
                else:
                    if len(infection_types) > 1:
                        is_non_cf_transmissible = True
            else:
                if "NA" in infection_types:
                    if len(infection_types) > 1:
                        is_non_cf_transmissible = True
                else:
                    if len(patients) > 1:
                        is_non_cf_transmissible = True
                
            print(node, infection_type, is_within_patient, is_cf_transmissible, is_non_cf_transmissible)
            st = tree.split("/")[3]
            out_type.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(node, st, infection_type, is_within_patient, is_cf_transmissible, is_non_cf_transmissible))
    out_type.close()
            



    








if __name__ == "__main__":
    parse()
