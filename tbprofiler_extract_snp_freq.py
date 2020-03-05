import json
import argparse
import os

def main(args):
    if args.itol:
        I = open("%s_%s.itol.txt" % (args.gene,args.mutation.replace(">","")),"w")
        I.write("""DATASET_SIMPLEBAR
SEPARATOR TAB
DATASET_LABEL\tmutation frequency
COLOR\t#ff0000

DATA
""")
    samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]
    for s in samples:
        j = json.load(open("%s/%s%s" % (args.dir,s,args.suffix)))
        mutation_not_found = True
        for var in j["variants"]:
            if var["gene_id"]==args.gene and var["change"]==args.mutation:
                print("%s\t%s" % (s,var["freq"]))
                if args.itol:
                    I.write("%s\t%s\n" % (s,var["freq"]))
                mutation_not_found = False
        if mutation_not_found:
            print("%s\t0.0" % s)
            if args.itol:
                I.write("%s\t%s\n" % (s,"0.0"))

    if args.itol:
        I.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('gene')
parser.add_argument('mutation')
parser.add_argument('--dir',type=str,default="results",help='Result directory')
parser.add_argument('--suffix',type=str,default=".results.json",help='Output file suffix')
parser.add_argument('--itol',action="store_true")
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
