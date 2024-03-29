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
import statistics

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
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    variant_freqs = []
    rprs = []
    bqrs = []
    mqrs = []
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            if (var["locus_tag"]==args.gene or var["gene"]==args.gene) and var["change"]==args.variant:
                variant_freqs.append(var["freq"])
                try:
                    rprs.append(float(var["variant_annotations"]["ReadPosRankSum"]))
                    bqrs.append(float(var["variant_annotations"]["BaseQRankSum"]))
                    mqrs.append(float(var["variant_annotations"]["MQRankSum"]))
                except:
                    pass

    if len(variant_freqs)>0:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (args.gene,args.variant,len(variant_freqs),statistics.median(variant_freqs),statistics.median(rprs),statistics.median(bqrs),statistics.median(mqrs)))
    else:
        print("%s\t%s\tNA\tNA\tNA\tNA\tNA" % (args.gene,args.variant))


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--gene',type=str,help='Gene',required=True)
parser.add_argument('--variant',type=str,help='Variant',required=True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
