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


def main(args):
    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    fields = ["sample_id","genome_pos","gene","change","freq","type","sublin","drtype","drugs"]
    with open(args.outfile,"w") as O:
        sep = "," if args.fmt == "csv" else "\t"
        writer = csv.DictWriter(O,fieldnames=fields,delimiter=sep)
        writer.writeheader()
        # Loop through the sample result files
        for s in tqdm(samples):
            # Data has the same structure as the .result.json files
            data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
            for var in data["dr_variants"] + data['other_variants']:
                var['sample_id'] = s
                var['sublin'] = data['sublin']
                var['drtype'] = data['drtype']
                if 'drugs' in var:
                    var['drugs'] = ",".join([a["drug"] for a in var["drugs"]])
                else:
                    var['drugs'] = ''
                writer.writerow({f:var[f] for f in fields})

# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--outfile',type=str,help='File with samples',required = True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--fmt',default="csv",choices=["csv","tsv"],type=str,help='Separator for fields')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
