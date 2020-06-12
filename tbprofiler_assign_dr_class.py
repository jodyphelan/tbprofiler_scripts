#! /usr/bin/env python

import sys
import csv
import json
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def main(args):
    bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix,args.db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    FLQ_set = set(["moxifloxacin","levofloxacin","ciprofloxacin","ofloxacin"])
    SLI_set = set(["amikacin","capreomycin","kanamycin"])

    OUT = open(args.out,"w")
    writer = csv.DictWriter(OUT, fieldnames = ["sample","dr-class"])
    writer.writeheader()

    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        resistant_drugs = set([d["drug"] for d in data["dr_variants"]])

        resistant_drugs = set([d["drug"] for d in data["dr_variants"]])
        rif = "rifampicin" in resistant_drugs
        inh = "isoniazid" in resistant_drugs
        flq = len(FLQ_set.intersection(resistant_drugs)) > 0
        sli = len(SLI_set.intersection(resistant_drugs)) > 0

        if len(resistant_drugs)==0:
            drtype = "Sensitive"
        elif (rif and not inh) or (inh and not rif):
            drtype = "Pre-MDR"
        elif (rif and inh) and (not flq and not sli):
            drtype = "MDR"
        elif (rif and inh) and ( (flq and not sli) or (sli and not flq) ):
            drtype = "Pre-XDR"
        elif (rif and inh) and (flq and sli):
            drtype = "XDR"
        else:
            drtype = "Other"

        writer.writerow({"sample":s, "dr-class":drtype})

    OUT.close()

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out', type=str, help='Name of CSV output', required=True)
parser.add_argument('--samples', type=str, help='File with samples')
parser.add_argument('--dir', default="results/", type=str, help='Directory containing results')
parser.add_argument('--db', default="tbdb", type=str, help='Database name')
parser.add_argument('--suffix', default=".results.json", type=str, help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
