#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict,Counter
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
    mapping = {
        "missense":"SNP",
        "non_coding":"SNP",
        "non_coding": "SNP",
        "stop_gained": "SNP",
        "start_lost": "SNP",
        "frameshift":"indel",
        "inframe_deletion":"indel",
        "inframe_insertion":"indel",
        "large_deletion":"large_deletion"
    }
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])

    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    resistance = defaultdict(lambda:defaultdict(list))
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        for var in data["dr_variants"]:
            resistance[var["drug"]][s].append(mapping.get(var["type"],"complex"))

    for drug in resistance:
        lines = []
        lines.append("DATASET_PIECHART")
        lines.append("SEPARATOR COMMA")
        lines.append("DATASET_LABEL,%s" % drug)
        lines.append("COLOR,#ff0000")
        lines.append("FIELD_COLORS,#ff0000,#00ff00,#0000ff,#ffffff")
        lines.append("FIELD_LABELS,snp,indel,large_deletion,no_variant")
        lines.append("MARGIN,5")
        # lines.append("MAXIMUM_SIZE,30")
        lines.append("BORDER_WIDTH,1")
        lines.append("BORDER_COLOR,#000000")
        lines.append("DATA")
        for s in samples:
            count = Counter(resistance[drug][s])
            lines.append("%s,-1,7,%s,%s" % (s, ",".join([str(count[d]) for d in ["SNP","indel","large_deletion"]]),"0" if sum(count.values())>0 else "1"))
        with open("%s.itol.conf.txt" % drug,"w") as O:
            O.write("\n".join(lines))




# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
