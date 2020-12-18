#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict, Counter
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
import statsmodels.api as sm
import numpy as np


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

def main(args):
    # Get a dictionary with the database file: {"ref": "/path/to/fasta" ... etc. }
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)

    # Get a dictionary mapping the locus_tags to drugs: {"Rv1484": ["isoniazid","pyrazinamide"], ... etc. }
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])

    locus_tag2genes = tbprofiler.rv2genes(conf["bed"])
    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    meta = {}
    for row in csv.DictReader(open(args.meta)):
        meta[row["run_accession"]] = row

    all_mutations = defaultdict(list)
    resistance_mutations = defaultdict(dict)

    results = []

    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"]+data["other_variants"]:
            if var["locus_tag"] not in locus_tag2genes: continue
            # if var["locus_tag"]=="Rv0676c": import pdb; pdb.set_trace()
            all_mutations[var["locus_tag"]].append(var["type"].replace("*",""))
            if "drugs" in var:
                for d in var["drugs"]:
                    resistance_mutations[d["drug"]][s] = var["type"].replace("*","")


    for drug in resistance_mutations:
        if drug not in list(meta.values())[0]: continue
        for vartype in set(resistance_mutations[drug].values()):
            tab = [
                [0, 0],
                [0, 0]
            ]

            sample_list = set(samples).intersection(set(list(meta)))
            tab[0][0] += len([s for s in sample_list if resistance_mutations[drug].get(s,"NA")==vartype and meta[s][drug]=="1"])
            tab[0][1] += len([s for s in sample_list if resistance_mutations[drug].get(s,"NA")!=vartype and meta[s][drug]=="1"])
            tab[1][0] += len([s for s in sample_list if resistance_mutations[drug].get(s,"NA")==vartype and meta[s][drug]=="0"])
            tab[1][1] += len([s for s in sample_list if resistance_mutations[drug].get(s,"NA")!=vartype and meta[s][drug]=="0"])

            t = sm.stats.Table2x2(np.asarray(tab))
            OR = t.oddsratio if (tab[0]!=[0.5,0.5] and (tab[0][0]!=0 and tab[1][0]!=0)) else "NA"
            OR_pval = t.oddsratio_pvalue() if (tab[0]!=[0.5,0.5] and (tab[0][0]!=0 and tab[1][0]!=0)) else "NA"
            pval = t.test_nominal_association().pvalue

            samples_with_vartype = sum([1 for s in resistance_mutations[drug] if resistance_mutations[drug][s]==vartype])
            samples_with_dr_var = len(resistance_mutations[drug])
            contribution_proportion =  samples_with_vartype/ samples_with_dr_var

            results.append({
                "drug":drug, "vartype":vartype,"oddsratio":OR, "oddsratio_pvalue": OR_pval, "pval": pval,
                "samples with vartype":samples_with_vartype,"samples_with_dr_var": samples_with_dr_var,
                "contribution_proportion": contribution_proportion
            })

        with open(args.out+".dr_variants.csv","w") as O:
            writer = csv.DictWriter(O,fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)

        with open(args.out+".all_variants.csv","w") as O:
            writer = csv.DictWriter(O,fieldnames=["gene","type","count"])
            writer.writeheader()
            for gene in all_mutations:
                counts = Counter(all_mutations[gene])
                for key in counts:
                    writer.writerow({"gene":locus_tag2genes[gene],"type":key,"count":counts[key]})


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out',type=str,help='File with samples')
parser.add_argument('--meta',type=str,help='File with samples')
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
