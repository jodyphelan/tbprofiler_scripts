#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
import copy

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

    # Get a dictionary mapping the locus_tags to drugs: {"Rv1484": ["isoniazid","ethionamide"], ... etc. }
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])

    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.in_dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        new_dr_variants = defaultdict(list)
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.in_dir,s,args.suffix))))
        for var in data["dr_variants"]:
            tmp = copy.deepcopy(var)
            del tmp["drug"]
            x = {"drug":var["drug"]}
            if "confidence" in tmp:
                del tmp["confidence"]
                x["confidence"] = var["confidence"]
            if "literature" in tmp:
                del tmp["literature"]
                x["literature"] = var["literature"]

            new_dr_variants[json.dumps(tmp)].append(x)

        data["dr_variants"] = []
        for x in new_dr_variants:
            new_var = json.loads(x)
            new_var["drugs"] = []
            for d in new_dr_variants[x]:
                new_var["drugs"].append(d)
            data["dr_variants"].append(new_var)
        json.dump(data,open("%s/%s%s" % (args.out_dir,s,args.suffix),"w"))


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--in-dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--out-dir',default="results_v3/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
