#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import re
from itertools import chain

def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

true_variants = {
    "por1": [("gyrA","p.Asp94Ala"),("rpoB","p.Ser450Leu"),("rrs","r.1401a>g"),("fabG1","c.-15C>T"),("inhA","p.Ile194Thr"),("pncA","p.Val125Gly"),("embA","c.-16C>T"),("embB","p.Met306Val"),("embB","p.Met423Thr"),("gid","p.Ala80Pro")],

}

def main(args):
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]
    sample_mutations = defaultdict(lambda:defaultdict(set))
    for s in tqdm(samples):
        biological_sample_name = ""
        re_obj = re.search("H37Rv[\d]_mq[\d]+_bq[\d]+",s)
        if re_obj:
            biological_sample_name = "H37Rv"
        re_obj = re.search("(por[\d]+).*_mq[\d]+_bq[\d]+",s)
        if re_obj:
            biological_sample_name = re_obj.group(1)


        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"]:
            if var["type"]=="missense":
                sample_mutations[biological_sample_name][s].add((var["gene"],var["change"]))


    for s in sample_mutations:
        union_mutations = set(list(chain(*[sample_mutations[s][run] for run in sample_mutations[s]])))
        for run in sample_mutations[s]:
            diff = union_mutations - set(sample_mutations[s][run])
            print("%s\t%s\t%s\t%s" % (s,run,len(union_mutations),len(sample_mutations[s][run])))





parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
