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
        samples = [x.replace(args.vcf_suffix,"") for x in os.listdir(args.vcf_dir) if x[-len(args.vcf_suffix):]==args.vcf_suffix]

    for l in open(conf["gff"]):
        row = l.strip().split()
        if len(row)<=2: continue
        if row[2]!="gene":continue
        if "Name=%s" % args.gene in l or "gene:%s" % args.gene in l:
            break

    start,end = int(row[3]),int(row[4])
    # Loop through the sample result files
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        if not os.path.isfile("%s/%s%s" % (args.dir,s,args.suffix)): continue
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        vars = json.dumps([d for d in data["dr_variants"]+data["other_variants"] if d["locus_tag"]==args.gene])
        print(vars)
        if "deletion" not in vars and "frameshift" not in vars and "inframe" not in vars and "stop" not in vars and "start" not in vars:

            revseq = "| revseq  -sequence /dev/stdin  -outseq /dev/stdout" if row[6]=="-" else ""
            pp.run_cmd("samtools faidx %s Chromosome:%s-%s | bcftools consensus %s/%s%s %s | sed 's/^>.*/>%s/' > %s.%s.fasta" % (conf["ref"],start,end,args.vcf_dir,s,args.vcf_suffix,revseq,s,s,args.gene),verbose=1)




# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--gene',type=str,help='File with samples',required=True)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--vcf-dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--vcf-suffix',default=".targets.vcf.gz",type=str,help='File suffix')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
