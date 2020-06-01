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
    drug2genes = defaultdict(set)
    gene2drugs = defaultdict(set)
    for l in open(conf["bed"]):
        row = l.strip().split()
        for drug in row[5].split(","):
            drug2genes[drug].add(row[4])
            gene2drugs[row[4]].add(drug)
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]
    drug_haplotypes = {drug:defaultdict(int) for drug in drug2genes}
    for s in tqdm(samples):
        sample_drug_haplotypes = {drug:set() for drug in drug2genes}
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"]:
            sample_drug_haplotypes[drug].add((var["gene"],var["change"]))
        for var in data["other_variants"]:
            if var["type"]=="synonymous": continue
            for drug in gene2drugs[var["gene"]]:
                sample_drug_haplotypes[drug].add((var["gene"],var["change"]))
        for drug in drug2genes:
            if len(sample_drug_haplotypes[drug])>0:
                drug_haplotypes[drug][tuple(sample_drug_haplotypes[drug])]+=1

    for drug in drug2genes:
        mutations = set(chain(*drug_haplotypes[drug]))
        for m1 in tqdm(mutations):
            for m2 in mutations:
                table = [[0,0],[0,0]]
                for haplotype in drug_haplotypes[drug]:
                    if m1 not in haplotype and m2 not in haplotype:
                        table[0][0] += drug_haplotypes[drug][haplotype]
                    if m1 in haplotype and m2 not in haplotype:
                        table[1][0] += drug_haplotypes[drug][haplotype]
                    if m1 not in haplotype and m2 in haplotype:
                        table[1][0] += drug_haplotypes[drug][haplotype]
                    if m1 in haplotype and m2 in haplotype:
                        table[1][1] += drug_haplotypes[drug][haplotype]
                print("%s_%s - %s_%s" % (m1[0],m1[1],m2[0],m2[1]), table)

        # for haplotype in drug_haplotypes[drug]:
        #     print("%s\t%s\t%s" % (
        #         drug,
        #         "__".join(["%s:%s" % (g,c) for g,c in haplotype]),
        #         drug_haplotypes[drug][haplotype]
        #     ))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
