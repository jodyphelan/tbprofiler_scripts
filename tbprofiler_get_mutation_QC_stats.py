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
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.results_dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    samples_with_mutation = []
    variant_position_set = set()
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.results_dir,s,args.suffix))))
        for var in data["dr_variants"] + data["other_variants"]:
            if (var["gene"]==args.gene or var["locus_tag"]==args.gene) and var["change"]==args.variant:
                samples_with_mutation.append(s)
                variant_position_set.add(var["genome_pos"])

    sys.stderr.write("\nFound %s samples with mutation\n" % len(samples_with_mutation))
    # samples_with_mutation = ["ERR2515541","ERR2510504","ERR2864225","SRR7341698"]
    if len(samples_with_mutation)==0:
        sys.stdout.write("%s\t%s\t%s\n" % (args.gene,args.variant,"Mutation_not_found"))
        quit()
    elif len(variant_position_set)>1:
        sys.stdout.write("%s\t%s\t%s\n" % (args.gene,args.variant,"Multiple_genome_pos"))
        quit()


    if len(variant_position_set)==1:
        variant_position = int(list(variant_position_set)[0])

    sys.stderr.write("\nGenome position is %s\n" % variant_position)
    sys.stderr.write("\nPerforming ReadPosRankSum test\n")
    # variant_position = 3841662
    params = vars(args)
    params["ref"] = conf["ref"]
    params["pos"] = variant_position
    params["tmp_vcf"] = pp.get_random_file(extension=".vcf.gz")
    read_pos_rank_sums = []
    for s in tqdm(samples_with_mutation):
        params["sample"] = s
        pp.run_cmd("tabix -f %(vcf_dir)s/%(sample)s.targets.csq.vcf.gz" % params,verbose=0)
        pp.run_cmd("bcftools view %(vcf_dir)s/%(sample)s.targets.csq.vcf.gz Chromosome:%(pos)s -Oz -o %(tmp_vcf)s" % params,verbose=0)
        pp.run_cmd("tabix -f %(tmp_vcf)s" % params,verbose=0)
        for l in pp.cmd_out("gatk VariantAnnotator -R %(ref)s -I %(bam_dir)s/%(sample)s%(bam_extension)s -V %(tmp_vcf)s -O /dev/stdout -A ReadPosRankSumTest -OVI false  | bcftools query -f '%%POS\\t%%ReadPosRankSum\\n'" % params,verbose=0):
            row = l.strip().split()
            if row[1]==".": continue
            if int(row[0])==variant_position:
                read_pos_rank_sums.append((s,float(row[1])))

    if len(read_pos_rank_sums)==0:
        sys.stdout.write("%s\t%s\t%s\n" % (args.gene,args.variant,"No_values_from_samples"))
    else:
        sys.stdout.write("%s\t%s\t%s\n" % (args.gene,args.variant,statistics.median([x[1] for x in read_pos_rank_sums])))
    pp.rm_files([params["tmp_vcf"]])

# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--gene',type=str,help='Directory containing results',required=True)
parser.add_argument('--variant',type=str,help='Directory containing vsfs',required=True)
parser.add_argument('--results-dir',type=str,help='Directory containing results',required=True)
parser.add_argument('--vcf-dir',type=str,help='Directory containing vsfs',required=True)
parser.add_argument('--bam-dir',type=str,help='Directory containing bams',required=True)
parser.add_argument('--bam-extension',default=".bqsr.cram",type=str,help='Extension to bams')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
