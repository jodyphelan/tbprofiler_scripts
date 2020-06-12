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

    blacklist = set([l.strip() for l in open(args.blacklist).readlines()]) if args.blacklist else []
    gene_set = set(args.gene.split(","))

    sys.stdout.write("sample,%s\n" % args.gene)
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        mutations = []
        for var in (data["dr_variants"] + data["other_variants"]):
            if var["gene"] in gene_set or var["locus_tag"] in gene_set:
                if not args.synonymous and var["type"]=="synonymous":
                    continue
                if var["change"] in blacklist:
                    continue
                mutations.append("%s_%s" % (var["gene"],var["change"]))

        sys.stdout.write("%s,%s\n" % (s,";".join(mutations)))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--gene',type=str,help='File with samples',required=True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--synonymous',action="store_true",help='File with samples')
parser.add_argument('--blacklist',type=str,help='Remove mutations from analysis')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
