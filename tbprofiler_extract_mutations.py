#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
from itertools import chain

def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf


def main(args):
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    drugs = set()
    rv2drugs = {}
    drugs2rv = {}
    for line in open(conf["bed"]):
        row = line.strip().split()
        rv2drugs[row[3]] = row[5].split(",")
        for drug in row[5].split(","):
            drugs.add(drug)
            drugs2rv[drug] = row[3]
    drugs = sorted(list(drugs))

    sys.stdout.write("sample,%s\n" % (",".join(["dr_mutations_%s,other_mutations_%s" % (d,d) for d in drugs])))
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        mutations = defaultdict(set)
        for var in data["dr_variants"]+data["other_variants"]:
            if var["type"]=="synonymous": continue
            if "drug" in var:
                mutations["dr_mutations_"+var["drug"]].add(var["gene"]+"_"+var["change"])
            else:
                tmp_drugs = rv2drugs[var["locus_tag"]]
                for drug in tmp_drugs:
                    mutations["other_mutations_"+drug].add(var["gene"]+"_"+var["change"])
        sys.stdout.write("%s,%s\n" % (s,",".join(["%s,%s" % (";".join(mutations["dr_mutations_"+d]),";".join(mutations["other_mutations_"+d])) for d in drugs])))



parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
