#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp


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

    for s in tqdm(samples):

        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        sys.stdout.write("%s\t%s\t%s\n" % (s,data["main_lin"],data["sublin"]))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".lineage.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
