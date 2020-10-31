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

    aa2genome_pos = defaultdict(list)
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))

        if args.dr_only:
            pool = data["dr_variants"]
        else:
            pool = data["dr_variants"] + data["other_variants"]
        for var in pool:
            aa2genome_pos[(var["gene"],var["change"])].append(var["genome_pos"])

    with open(args.out,"w") as O:
        for var in aa2genome_pos:
            for pos in aa2genome_pos[var]:
                O.write("%s\t%s\t%s\n" % (var[0],var[1],pos))

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out',type=str,help='Output file',required = True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.add_argument('--dr-only',action="store_true",help="Only extract drug resistance associated variants")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
