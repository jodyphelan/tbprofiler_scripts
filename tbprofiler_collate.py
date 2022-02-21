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
    files = {'gff':'.gff','ref':'.fasta','ann':'.ann.txt','barcode':'.barcode.bed','bed':'.bed','json_db':'.dr.json','version':'.version.json'}
    conf = {}
    for key in files:
        sys.stderr.write('Using %s file: %s\n' % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf


def main(args):
    # Get a dictionary with the database file: {'ref': '/path/to/fasta' ... etc. }
    conf = get_conf_dict(sys.base_prefix + '/share/tbprofiler/%s' % args.db)

    # Get a dictionary mapping the locus_tags to drugs: {'Rv1484': ['isoniazid','ethionamide'], ... etc. }
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf['bed'])
    
    # Get a dictionary mapping the drug to genes: {'rifampicin': ['rpoB', 'rpoC'], 'clofazimine': ['mmpR5', 'pepQ'], ... etc. }
    drug2genes = tbprofiler.get_drugs2gene(conf['bed'])
    

    # If a list of samples is supplied through the args object, store it in a list else get the list from looking in the results direcotry
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,'') for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    # Loop through the sample result files
    drugs = [
        'rifampicin','isoniazid','ethambutol','pyrazinamide','streptomycin','amikacin',
        'kanamycin','capreomycin','fluoroquinolones','ethionamide','cycloserine',
        'para-aminosalicylic_acid','clofazimine','bedaquiline','delamanid'
    ]

    # Set up a list which will contain our output file rows
    rows = []

    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck(f'{args.dir}/{s}{args.suffix}')))

        # The data is organised per variant in data['dr_variants']. We need to 
        # transform this into a structure which is arranged by drug instead.
        # We do this by:
        # 1. Setting up a dictionary (drug_variants) where the values are lists
        # 2. Loop through all the variants and append the gene/change/freq
        #    to the list for each drug
        #
        # The structure will look like {'isoniazid':['katG_p.Ser315Thr_0.95','fabG1_-15T>C_1.00']}
        drug_variants = defaultdict(list)
        for var in data['dr_variants']:
            for d in var['drugs']:
                drug_variants[d['drug']].append(f'{var["gene"]}_{var["change"]}_{round(var["freq"],2)}')

        # Create a lookup dictionary containing all the genes for which we have missing coverage
        # E.g. {'rpoB':'0.2'}
        gene_coverage = {d['gene']:str(d['fraction']) for d in data['qc']['gene_coverage'] if d['fraction']>0}

        # Set up our row for the final output file with column names being the keys.
        row = {
            'sample': s,
            'main_lineage': data['main_lin'],
            'sublineage': data['sublin'],
            'drtype': data['drtype']
        }
        
        # For each drug add a column with the variants and another containing a value if there is missin coverage
        for drug in drugs:
            row[f'{drug}_variants'] = ", ".join(drug_variants[drug])
            row[f'{drug}_gene_cov'] = ", ".join([gene_coverage[gene] for gene in drug2genes[drug] if gene in gene_coverage])
        rows.append(row)


    # Write the output file
    with open(args.outfile,'w') as O:
        writer = csv.DictWriter(O,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--outfile',type=str,help='Output file')
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',default='results/',type=str,help='Directory containing results')
parser.add_argument('--db',default='tbdb',type=str,help='Database name')
parser.add_argument('--suffix',default='.results.json',type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
