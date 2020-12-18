#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
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
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    mutations = defaultdict(lambda: defaultdict(int))
    gene_set = set(args.genes.split(",")) if args.genes else None
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var1 in (data["dr_variants"] + data["other_variants"]):
            if args.genes and (var1["gene"] not in gene_set and var1["locus_tag"] not in gene_set):
                continue
            if not args.synonymous and var1["type"]=="synonymous":
                continue
            for var2 in (data["dr_variants"] + data["other_variants"]):
                if args.genes and (var2["gene"] not in gene_set and var2["locus_tag"] not in gene_set):
                    continue
                if not args.synonymous and var2["type"]=="synonymous":
                    continue
                mutations[var1["gene"]+"_"+var1["change"]][var2["gene"]+"_"+var2["change"]] += 1

    for mut1 in tqdm(mutations):
        for mut2 in mutations[mut1]:
            total_mut1 = mutations[mut1][mut1]
            total_mut2 = mutations[mut2][mut2]


            t = [
                [0,0],
                [0,0]
            ]
            t[0][0] = mutations[mut1][mut2]
            t[1][0] = total_mut1 - t[0][0]
            t[0][1] = total_mut2 - t[0][0]
            t[1][1] = len(samples) - t[0][1] - t[1][0] - t[0][0]
            t2 = sm.stats.Table2x2(np.asarray(t))
            OR = t2.oddsratio if t[0]!=[0.5,0.5] else "NA"
            OR_pval = t2.oddsratio_pvalue() if t[0]!=[0.5,0.5] else "NA"
            pval = t2.test_nominal_association().pvalue
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mut1,mut2,total_mut1,total_mut2,mutations[mut1][mut2],OR,OR_pval,pval),t)







parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--genes',type=str,help='Genes to include')
parser.add_argument('--synonymous',action="store_true",help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
