#! /usr/bin/env python
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import pathogenprofiler as pp
import csv

try:
    sys.base_prefix
except:
    sys.base_prefix = getattr(sys, 'base_prefix', getattr(sys, 'real_prefix', sys.prefix))


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf



def main(args):
    if args.drugs:
        args.drugs = [x.lower() for x in args.drugs.split(",")]
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(".results.json","") for x in os.listdir(args.dir) if x[-13:]==".results.json"]

    variants = defaultdict(set)
    variant_drug_associations = defaultdict(set)
    drtypes = defaultdict(set)
    genome_pos2changes = defaultdict(set)
    for s in tqdm(samples):
        tmp = json.load(open(f"{args.dir}/{s}.results.json"))
        drtypes[tmp["drtype"]].add(s)
        for var in tmp["dr_variants"]+tmp["other_variants"]:
            variants[var["genome_pos"]].add(s)
            genome_pos2changes[var["genome_pos"]].add(var["change"])

    total_sample_n = len(samples)
    OUT = open(args.outfile,"w")
    OUT.write("pos\tchange\ttotal_num\ttotal_af\t%s\t%s\n" % (
        "\t".join(["%s_num" % d for d in  drtypes]),
        "\t".join(["%s_af" % d for d in  drtypes])
    ))
    for pos in tqdm(sorted(variants)):
        total_num = len(variants[pos])
        af = total_num/total_sample_n
        sub_af = {}
        sub_num = {}
        for drtype in drtypes:
            mut = len([s for s in drtypes[drtype] if s in variants[pos]])
            no_mut = len([s for s in drtypes[drtype] if s not in variants[pos]])
            sub_num[drtype] = mut
            sub_af[drtype] = mut/(mut+no_mut) if (mut+no_mut)>0 else "NA"
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
            pos,
            ",".join(genome_pos2changes[pos]),
            total_num,
            af,
            "\t".join([str(sub_num[d]) for d in drtypes]),
            "\t".join([str(sub_af[d]) for d in drtypes])
        ))
    OUT.close()
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--outfile',type=str,help='NGS Platform',required=True)
parser.add_argument('--samples',type=str,help='NGS Platform')
parser.add_argument('--dir',default="results/",type=str,help='NGS Platform')
parser.add_argument('--db',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--drugs',type=str,help='NGS Platform')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
