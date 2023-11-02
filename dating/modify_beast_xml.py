import xml.etree.ElementTree as ET
import argparse
from Bio import SeqIO
import random 

TEMPLATE_RUN_ID ="313_wo_hyp_snp-sites" 
def parse():
    parser = argparse.ArgumentParser("read in BEAST xml and modify")
    parser.add_argument("beast_xml", help="input xml")
    parser.add_argument("run_id", help="id for this BEAST run")
    parser.add_argument("--dates", help="isolate sampling dates")
    parser.add_argument("--ucld_mean_uniform", help="set clock rate and lower and upper bound for uniform ucld mean prior (two order of magnitude on either side)")
    parser.add_argument("--alignment", help="alignment")
    parser.add_argument("--randomise_dates", help="randomise date trait", action = "store_true")
    parser.add_argument("--ucld_mean", help="set clock.rate and ucld mean log normal prior")
    parser.add_argument("--ucld_mean_std", help="set ucld log normal standard deviation prior")
    parser.add_argument("--chain_length", help="set chain length")
    parser.add_argument("--locations", help="sampling locations")
    parser.add_argument("xml_out", help="output xml")
    args = parser.parse_args()
    process(**vars(args))

def process(beast_xml, run_id, dates, alignment, randomise_dates, ucld_mean, ucld_mean_std, ucld_mean_uniform, chain_length, locations, xml_out):
    with open(beast_xml) as f:
        xml_string = f.read()
    xml_string = xml_string.replace(TEMPLATE_RUN_ID, run_id)
    if locations:
        xml_string = xml_string.replace("location", "location_%s" % run_id)
    xml_tree = ET.ElementTree(ET.fromstring(xml_string))
    xml_root = xml_tree.getroot()
    if alignment:
        data_node = xml_root.findall('data')[0]
        for seq in data_node.findall('sequence'):
            data_node.remove(seq)
        records = SeqIO.parse(alignment, "fasta")
        #TODO update rate categories = 2N - 2
        for record in records:
            seq_node = ET.SubElement(data_node, "sequence")
            seq_node.attrib = {"id": "seq_" + record.id, "spec": "Sequence", "taxon": record.id, "totalcount": "4", "value": record.seq}

    if ucld_mean_uniform:
        mean_rate_prior = xml_root.findall('''.//*[@id='MeanRatePrior.c:%s']''' % run_id)[0]
        if "," in ucld_mean_uniform:
            lower, upper = ucld_mean_uniform.split(',')
            mean_rate_prior[0].attrib["lower"] = lower 
            mean_rate_prior[0].attrib["upper"] = upper 
        else:
            mean_rate_prior[0].attrib["lower"] = str(ucld_mean_uniform / 100)
            mean_rate_prior[0].attrib["upper"] = str(ucld_mean_uniform * 100)
    
    if ucld_mean or ucld_mean_std:
        mean_rate_prior = xml_root.findall('''.//*[@id='MeanRatePrior.c:%s']''' % run_id)[0]
        #impo/t pdb
        #pdb.set_trace()
        for param in mean_rate_prior[0].findall('parameter'):
            if ucld_mean is not None and param.attrib["name"] == "M":
                param.text = str(ucld_mean)
            if ucld_mean_std is not None and param.attrib["name"] == "S":
                param.text = str(ucld_mean_std)

    if ucld_mean:
        ucld_mean_param = xml_root.findall('''.//*[@id='ucldMean.c:%s']''' % run_id)[0]
        ucld_mean_param.text = str(ucld_mean)
    if ucld_mean_uniform:
        ucld_mean_param = xml_root.findall('''.//*[@id='ucldMean.c:%s']''' % run_id)[0]
        ucld_mean_param.text = str((float(lower) + float(upper)) / 2)
    #import pdb
    #pdb.set_trace()
    if not chain_length is None:
        xml_root.findall('run')[0].attrib['chainLength'] = chain_length
    if dates:
        #parse dates
        isol2date = {}
        with open(dates, 'r') as f:
            f.readline()
            for l in f:
                sample, date = l.strip().split("\t")
                isol2date[sample] = date
        dates_node = xml_root.findall('''.//*[@id='dateTrait.t:%s']''' % run_id)[0]
        data_node = xml_root.findall('data')[0]
        isol_dates = []
        for seq in data_node.findall('sequence'):
            taxon = seq.attrib["taxon"]
            date = isol2date[taxon] 
            isol_dates.append("%s=%s" % (taxon, date))
        dates_node.attrib["value"] = ",".join(isol_dates)
    
    if locations:
        #parse locations 
        isol2location = {}
        all_locations = []
        with open(locations, 'r') as f:
            f.readline()
            for l in f:
                sample, location = l.strip().split("\t")
                isol2location[sample] = location 
        #find locations node
        locations_node = xml_root.findall(""".//*[@id='traitSet.location_%s']""" % run_id)[0]
        data_node = xml_root.findall('data')[0]
        isol_locations = []
        for seq in data_node.findall('sequence'):
            taxon = seq.attrib["taxon"]
            location = isol2location[taxon] 
            all_locations.append(location)
            isol_locations.append("%s=%s" % (taxon, location))
        locations_node.text = ",".join(isol_locations)
        #determine number of unique continents
        uq_locations = sorted(list(set(all_locations)))
        #update codemap
        codemap = xml_root.findall(""".//*[@id='traitDataType.location_%s']""" % run_id)[0]
        codemap_str =  ["{}={}".format(i, j) for i, j in zip(uq_locations, range(len(uq_locations)))]
        codemap_str = ",".join(codemap_str)
        codemap_str = codemap_str + ",? =" + " ".join([str(i) for i in range(len(uq_locations))]) + " "
        codemap.attrib['codeMap'] = codemap_str
        codemap.attrib['states'] = str(len(uq_locations))
        #update trait frequencies
        trait_freqs = xml_root.findall(""".//*[@id='traitfrequencies.s:location_%s']""" % run_id)[0]
        trait_freqs.attrib['dimension'] = str(len(uq_locations))
        trait_freqs.text = str(1/len(uq_locations))
        #offset Poisson
        poisson = xml_root.findall(""".//*[@id='Poisson.2']""")[0]
        poisson.attrib['offset'] = str(len(uq_locations) - 1)
        #update geo rates
        geo_rates = xml_root.findall(""".//*[@id='relativeGeoRates.s:location_%s']""" % run_id)[0]
        geo_rates.attrib['dimension'] = str(len(uq_locations) * (len(uq_locations) - 1)// 2)
        rate_indicator = xml_root.findall(""".//*[@id='rateIndicator.s:location_%s']""" % run_id)[0]
        rate_indicator.attrib['dimension'] = str(len(uq_locations) * (len(uq_locations) - 1)// 2)



    if randomise_dates:
        dates_node = xml_root.findall('''.//*[@id='dateTrait.t:%s']''' % run_id)[0]
        date_dict = dict([i.split("=") for i in dates_node.attrib["value"].split(",")])
        date_keys = [i for i in date_dict.keys()]
        random.shuffle(date_keys)
        dates_shuffled = ["%s=%s" % (key, value) for key, value in zip(date_keys, date_dict.values())]
        dates_node.attrib["value"] = ",".join(dates_shuffled)
    xml_tree.write(xml_out)

if __name__ == "__main__":
    parse()
