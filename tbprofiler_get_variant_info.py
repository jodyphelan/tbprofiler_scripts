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
import re

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
    mutations = defaultdict(list)
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            if var["gene"]!="rrs" and var["gene"]!="rrl" and re.search("r\.[0-9]+",var["change"]):
                continue
            if "nucleotide_change" in var and (re.search("p.[A-Za-z]+",var["change"]) or re.search("c.[0-9]+",var["change"]) or re.search("c.\-[0-9]+",var["change"])):
                pos =  ",".join([re.search("([0-9]+)([ACGT]+)>([ACGT]+)",x).group(1) for x in var["nucleotide_change"].split("+")])
                ref =  ",".join([re.search("([0-9]+)([ACGT]+)>([ACGT]+)",x).group(2) for x in var["nucleotide_change"].split("+")])
                alt =  ",".join([re.search("([0-9]+)([ACGT]+)>([ACGT]+)",x).group(3) for x in var["nucleotide_change"].split("+")])
                # if var["change"]=="p.Gly168Ser":
                    # import pdb; pdb.set_trace()
            elif var["type"]=="non_coding" and re.search("c.\-[0-9]+",var["change"]):
                re_obj = re.search("c.\-[0-9]+([ACGT]+)>([ACGT]+)",var["change"])
                pos = str(var["genome_pos"])
                ref = re_obj.group(1)
                alt = re_obj.group(2)
            elif var["type"]=="non_coding" and re.search("r.[0-9]+",var["change"]):
                re_obj = re.search("[0-9]+([ACGT]+)>([ACGT]+)", var["_internal_change"])
                pos = str(var["genome_pos"])
                ref = re_obj.group(1)
                alt = re_obj.group(2)
            elif var["type"]=="large_deletion":
                continue
            elif var["type"].replace("*","")=="synonymous":
                continue
            elif var["type"].replace("*","")=="frameshift&start_lost":
                continue
            elif var["type"].replace("*","")=="missense&inframe_altering":
                continue
            elif var["type"].replace("*","")=="stop_lost":
                continue
            elif var["type"].replace("*","")=="stop_retained":
                continue
            else:
                quit(var)
            # if var["change"]=="p.Ser450Leu":
                # import pdb; pdb.set_trace()
            mutations[(pos,ref,alt)].append(json.dumps({
                "genome_pos": pos,
                "type":var["type"].replace("*",""),
                "locus_tag":var["locus_tag"],
                "gene":var["gene"],
                "_internal_change":var["_internal_change"],
                "change":var["change"]
            }))

    for key in sorted(mutations,key=lambda x:int(x[0].split(",")[0])):
        for x in set(mutations[key]):
            var = json.loads(x)
            print("Chromosome\t%s\t%s\t%s\t%s|%s|%s|%s|%s|%s|%s" % (var["genome_pos"],key[1],key[2],var["type"],var["locus_tag"],var["gene"],"NA","NA",var["_internal_change"],var["change"]))


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
