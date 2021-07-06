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

def camel_case(c):
    d = list(c.replace("_"," ").replace("-", " ").title().replace(" ",""))
    d[0] = d[0].lower()
    d = "".join(d)
    return d

def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","ann":".ann.txt","barcode":".barcode.bed","bed":".bed","json_db":".dr.json","version":".version.json"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf

def uniq_dict_list(l):
    return [json.loads(j) for j in set([json.dumps(d) for d in l])]

def standardise_types(dictlist):
    new_dict_list = []
    id2type = defaultdict(set)
    for row in dictlist:
        id2type[row["id"]].add(row["type"])
    for row in dictlist:
        if id2type[row["id"]]==set(["frameshift","inframe_insertion"]):
            if row["type"]=="inframe_insertion": continue
        elif id2type[row["id"]]==set(["frameshift","inframe_deletion"]):
            if row["type"]=="inframe_deletion": continue
        new_dict_list.append(row)
    return new_dict_list

def main(args):
    conf = get_conf_dict(sys.base_prefix + "/share/tbprofiler/%s" % args.db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(conf["bed"])
    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    if args.meta:
        meta = {}
        for row in csv.DictReader(open(args.meta)):
            meta[row["wgs_id"]] = row



    sample_nodes = []
    tmp_variant_nodes = []
    sample_variant_edges = []
    drugs = set()
    lineage_nodes = set()
    tmp_variant_drug_edges = []
    sample_lineage_edges = []
    spoligotype_nodes = set()
    sample_spoligotype_edges = []
    if args.spoligotypes:
        spoligotypes = {}
        for row in csv.DictReader(open(args.spoligotypes)):
            spoligotypes[row["sample"]] = row["spoligotype"]
            sample_spoligotype_edges.append({"id":row["sample"],"spoligotype":row["spoligotype"]})
            spoligotype_nodes.add(row["spoligotype"])

    # Loop through the sample result files
    for s in tqdm(samples):
        # Data has the same structure as the .result.json files
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        sample_node = {
            "id":s,
            "drtype":data["drtype"],
            "lineage":data["sublin"],
            "lineageInfo": json.dumps(data["lineage"]),
            "qc": json.dumps({"pct_reads_mapped":data["qc"]["pct_reads_mapped"],"pct_reads_mapped":data["qc"]["pct_reads_mapped"],"gene_coverage":[]}),
            "pipeline": json.dumps(data["pipeline"]),
            "tbprofilerVersion": json.dumps(data["tbprofiler_version"]),
            "dbVersion": json.dumps(data["db_version"]),
        }
        lineage_nodes.add(data["sublin"])
        sample_lineage_edges.append({"sampleId":s,"lineage":data["sublin"]})
        if args.meta:
            for c in list(meta.values())[0]:
                if c=="id": continue
                d = camel_case(c)
                sample_node[d] = meta[s][c] if s in meta else "NA"

        if args.spoligotypes:
            sample_node["spoligotype"] = spoligotypes.get(s,"NA")

        sample_nodes.append(sample_node)
        for var in data["dr_variants"] + data["other_variants"]:
            # if var["type"]=="synonymous": continue
            variant_id = "%s_%s" % (var["locus_tag"],var["change"])
            sample_variant_edges.append(
                {
                    "sampleId":s,
                    "variantId":variant_id,
                    "freq":var["freq"],
                    "genome_pos": var["genome_pos"],
                    "nucleotideChange": var.get("nucleotide_change","NA"),
                    "internalChange": var.get("_internal_change","NA")
                }
            )
            tmp_variant_nodes.append(
                {
                    "id": variant_id,
                    "type": var["type"].replace("*",""),
                    "change": var["change"],
                    "gene": var["gene"],
                    "locus_tag": var["locus_tag"],
                }
            )
            if "drugs" in var:
                for d in var["drugs"]:
                    drugs.add(d["drug"])
                    tmp_variant_drug_edges.append(
                        {
                            "variantId": variant_id,
                            "drug": d["drug"]
                        }
                    )


    drug_nodes = [{"id":d} for d in drugs]
    lineage_nodes = [{"id":d} for d in lineage_nodes]
    spoligotype_nodes = [{"id":d} for d in spoligotype_nodes]
    variant_drug_edges = uniq_dict_list(tmp_variant_drug_edges)
    variant_nodes = uniq_dict_list(standardise_types(tmp_variant_nodes))

    def batch(iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]

    with open("sample_nodes.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(sample_nodes[0]))
        writer.writeheader()
        writer.writerows(sample_nodes)

    with open("variant_nodes.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(variant_nodes[0]))
        writer.writeheader()
        writer.writerows(variant_nodes)

    for i,x in enumerate(batch(list(range(len(sample_variant_edges))),10000)):
        with open("sample_variant_edges.%s.csv" % i,"w") as O:
            writer = csv.DictWriter(O,fieldnames = list(sample_variant_edges[0]))
            writer.writeheader()
            for j in x:
                writer.writerow(sample_variant_edges[j])

    with open("drug_nodes.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(drug_nodes[0]))
        writer.writeheader()
        writer.writerows(drug_nodes)

    with open("variant_drug_edges.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(variant_drug_edges[0]))
        writer.writeheader()
        writer.writerows(variant_drug_edges)

    with open("lineage_nodes.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(lineage_nodes[0]))
        writer.writeheader()
        writer.writerows(lineage_nodes)

    with open("sample_lineage_edges.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(sample_lineage_edges[0]))
        writer.writeheader()
        writer.writerows(sample_lineage_edges)

    with open("spoligotype_nodes.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(spoligotype_nodes[0]))
        writer.writeheader()
        writer.writerows(spoligotype_nodes)

    with open("sample_spoligotype_edges.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames = list(sample_spoligotype_edges[0]))
        writer.writeheader()
        writer.writerows(sample_spoligotype_edges)

    with open("Cypher_commands.txt" ,"w") as O:
        if args.index:
            O.write("CREATE INDEX FOR (n:Sample) ON (n.id);\n")
            O.write("CREATE INDEX FOR (n:SRA) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///sample_nodes.csv' AS csvLine\n")
        O.write("CREATE (s:Sample:SRA {%s});\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in sample_nodes[0]]))
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Country) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///sample_nodes.csv' AS csvLine\n")
        O.write("WITH csvLine WHERE NOT csvLine.countryCode IS null\n")
        O.write("MERGE (s:Sample {id:csvLine.id})\n")
        O.write("MERGE (c:Country {id:csvLine.countryCode})\n")
        O.write("CREATE (s) -[:COLLECTED_IN]-> (c);\n")
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Variant) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///variant_nodes.csv' AS csvLine\n")
        O.write("CREATE (v:Variant {%s});\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in variant_nodes[0]]))
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Gene) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///variant_nodes.csv' AS csvLine\n")
        O.write("MERGE (v:Variant {id:csvLine.id})\n")
        O.write("MERGE (g:Gene {id:csvLine.locus_tag, locusTag:csvLine.locus_tag, name:csvLine.gene})\n")
        O.write("CREATE (v) -[:IN_GENE]-> (g);\n")
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Drug) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///drug_nodes.csv' AS csvLine\n")
        O.write("CREATE (d:Drug {%s});\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in drug_nodes[0]]))
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Lineage) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///lineage_nodes.csv' AS csvLine\n")
        O.write("CREATE (:Lineage {%s});\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in lineage_nodes[0]]))
        O.write("\n")

        for i,x in enumerate(batch(list(range(len(sample_variant_edges))),10000)):
            O.write("LOAD CSV WITH HEADERS FROM 'file:///sample_variant_edges.%s.csv' AS csvLine\n" % i)
            O.write("MATCH (s:Sample {id: csvLine.sampleId}),(v:Variant {id:csvLine.variantId})\n")
            O.write("CREATE (s) -[:CONTAINS {%s}]-> (v);\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in sample_variant_edges[0]]))
            O.write("\n")

        O.write("LOAD CSV WITH HEADERS FROM 'file:///variant_drug_edges.csv' AS csvLine\n")
        O.write("MATCH (v:Variant {id: csvLine.variantId}),(d:Drug {id:csvLine.drug})\n")
        O.write("CREATE (v) -[:CONFERS_RESISTANCE]-> (d);\n")
        O.write("\n")

        O.write("LOAD CSV WITH HEADERS FROM 'file:///sample_lineage_edges.csv' AS csvLine\n")
        O.write("MATCH (s:Sample {id: csvLine.sampleId}),(l:Lineage {id:csvLine.lineage})\n")
        O.write("CREATE (s) -[:LINEAGE]-> (l);\n")
        O.write("\n")

        if args.index:
            O.write("CREATE INDEX FOR (n:Spoligotype) ON (n.id);\n")
        O.write("LOAD CSV WITH HEADERS FROM 'file:///spoligotype_nodes.csv' AS csvLine\n")
        O.write("CREATE (:Spoligotype {%s});\n" % ", ".join(["%s: csvLine.%s" % (d,d) for d in spoligotype_nodes[0]]))
        O.write("\n")

        O.write("LOAD CSV WITH HEADERS FROM 'file:///sample_spoligotype_edges.csv' AS csvLine\n")
        O.write("MATCH (s:Sample {id: csvLine.id}),(l:Spoligotype {id:csvLine.spoligotype})\n")
        O.write("CREATE (s) -[:SPOLIGOTYPE]-> (l);\n")
        O.write("\n")


        # CREATE CONSTRAINT ON (m:variant) ASSERT variant.id IS UNIQUE
# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--dir',type=str,help='Directory containing results',required=True)
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.add_argument('--meta',type=str,help='Meta file')
parser.add_argument('--spoligotypes',type=str,help='Meta file',required=True)
parser.add_argument('--index',action="store_true",help='Write index commands')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
