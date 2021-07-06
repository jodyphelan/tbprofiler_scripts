#! /usr/bin/env python3
import sys
import pathogenprofiler as pp
import argparse
import json
import tbprofiler as tbp
import os
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

def main_profile(args):
    #### Setup conf dictionary ###
    if args.db=="tbdb" and not args.external_db and pp.nofile(sys.base_prefix+"/share/tbprofiler/tbdb.fasta"):
        pp.log("Can't find the tbdb file at %s. Please run 'tb-profiler update_tbdb' to load the default library or specify another using the '--external_db' flag" % sys.base_prefix,ext=True)
    if args.external_db:
        conf = get_conf_dict(args.external_db)
    else:
        conf = get_conf_dict(sys.base_prefix+"/share/tbprofiler/%s" % args.db)

    ### Create folders for results if they don't exist ###
    if pp.nofolder(args.dir):
        os.mkdir(args.dir)

    for x in ["bam","vcf","results"]:
        if pp.nofolder(args.dir+"/"+x):
            os.mkdir(args.dir+"/"+x)

    ### Set up platform dependant parameters ###
    if args.platform=="nanopore":
        args.mapper = "minimap2"
        args.caller = "bcftools"
        args.no_trim=True
        run_delly = False
    else:
        if args.no_delly:
            run_delly = False
        else:
            run_delly = True

    ### Setup prefix for files ###
    files_prefix = args.dir+"/"+args.prefix

    ### Create bam file if fastq has been supplied ###
    if args.bam==None:
        if args.read1 and args.read2 and args.no_trim:
            # Paired + no trimming
            fastq_obj = pp.fastq(args.read1,args.read2)
        elif args.read1 and args.read2 and not args.no_trim:
            # Paired + trimming
            untrimmed_fastq_obj = pp.fastq(args.read1,args.read2)
            fastq_obj = untrimmed_fastq_obj.trim(files_prefix,threads=args.threads)
        elif args.read1 and not args.read2 and args.no_trim:
            # Unpaired + trimming
            fastq_obj = pp.fastq(args.read1,args.read2)
        elif args.read1 and not args.read2 and not args.no_trim:
            # Unpaired + trimming
            untrimmed_fastq_obj = pp.fastq(args.read1)
            fastq_obj = untrimmed_fastq_obj.trim(files_prefix,threads=args.threads)
        else:
            exit("\nPlease provide a bam file or a fastq file(s)...Exiting!\n")
        bam_obj = fastq_obj.map_to_ref(
            ref_file=conf["ref"], prefix=files_prefix,sample_name=args.prefix,
            aligner=args.mapper, platform=args.platform, threads=args.threads
        )
        bam_file = bam_obj.bam_file
    else:
        bam_file = args.bam

    print(args.delly_bcf_file)
    run_coverage = False if args.no_coverage else True
    ### Run profiling module from pathogen-profiler ###
    results = pp.bam_profiler(
        conf=conf, bam_file=bam_file, prefix=files_prefix, platform=args.platform,
        caller=args.caller, threads=args.threads, no_flagstat=args.no_flagstat,
        run_delly = run_delly, calling_params=args.calling_params,
        coverage_fraction_threshold=args.coverage_fraction_threshold,
        missing_cov_threshold=args.missing_cov_threshold,
        delly_bcf_file=args.delly_bcf_file
    )
    json.dump(results,open(args.prefix+".tmp_results.json","w"))
    ### Reformat the results to TB-Profiler style ###
    results = tbp.reformat(results, conf, reporting_af=args.reporting_af)
    results["id"] = args.prefix
    results["tbprofiler_version"] = tbp._VERSION
    results["pipeline"] = {"mapper":args.mapper if not args.bam else "N/A","variant_caller":args.caller}

    json_output = args.dir+"/results/"+args.prefix+".results.json"
    tex_output = args.dir+"/results/"+args.prefix+".results.tex"
    text_output = args.dir+"/results/"+args.prefix+".results.txt"
    csv_output = args.dir+"/results/"+args.prefix+".results.csv"

    json.dump(results,open(json_output,"w"))
    extra_columns = [x.lower() for x in args.add_columns.split(",")] if args.add_columns else []
    if args.pdf:
        tbp.write_tex(results,conf,tex_output,extra_columns)
        pp.run_cmd("pdflatex %s"%tex_output,verbose=1)
        pp.rm_files([tex_output, args.dir+"/"+args.prefix+".results.aux",args.dir+"/"+args.prefix+".results.log"])
    if args.txt:
        tbp.write_text(results,conf,text_output,extra_columns,reporting_af=args.reporting_af)
    if args.csv:
        tbp.write_csv(results,conf,csv_output,extra_columns)

    ### Move files to respective directories ###
    if not args.bam:
        pp.run_cmd("mv %(dir)s/%(prefix)s.bam* %(dir)s/bam/" % vars(args))
        if not args.no_trim:
            pp.run_cmd("rm -f %s" % " ".join(fastq_obj.files))
    pp.run_cmd("mv -f %(dir)s/%(prefix)s*.vcf.gz* %(dir)s/vcf/" % vars(args))
    if run_delly and results["delly"]=="success" and not args.delly_bcf_file:
        pp.run_cmd("mv -f %(dir)s/%(prefix)s.delly.bcf* %(dir)s/vcf/" % vars(args))

    ### Add meta data to results
    if args.meta:
        for row in csv.DictReader(open(args.meta)):
            if row["id"]==results["id"]:
                for col in row:
                    results["meta_"+col] = row[col]
    pp.log("Profiling finished sucessfully!")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version="TBProfiler version %s" % tbp._VERSION)
    subparsers = parser.add_subparsers(help="Task to perform")

    parser_sub = subparsers.add_parser('profile', help='Run whole profiling pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_sub.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
    parser_sub.add_argument('--read1','-1',help='First read file')
    parser_sub.add_argument('--read2','-2',help='Second read file')
    parser_sub.add_argument('--bam','-a',help='BAM file. Make sure it has been generated using the H37Rv genome (GCA_000195955.2)')
    parser_sub.add_argument('--prefix','-p',default="tbprofiler",help='Sample prefix for all results generated')
    parser_sub.add_argument('--no_trim',action="store_true",help="Don't trim files using trimmomatic")
    parser_sub.add_argument('--db',default='tbdb',help='Mutation panel name')
    parser_sub.add_argument('--external_db',type=str,help='Path to db files prefix (overrides "--db" parameter)')
    parser_sub.add_argument('--mapper',default="bwa", choices=["bwa","minimap2","bowtie2"],help="Mapping tool to use. If you are using nanopore data it will default to minimap2",type=str)
    parser_sub.add_argument('--caller',default="bcftools", choices=["bcftools","gatk","freebayes"],help="Variant calling tool to use.",type=str)
    parser_sub.add_argument('--calling_params',type=str,help='Override default parameters for variant calling')
    parser_sub.add_argument('--min_depth',default=10,type=int,help='Minimum depth required to call variant. Bases with depth below this cutoff will be marked as missing')
    parser_sub.add_argument('--af',default=0.1,type=float,help='Minimum allele frequency to call variants')
    parser_sub.add_argument('--reporting_af',default=0.1,type=float,help='Minimum allele frequency to use variants for prediction')
    parser_sub.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
    parser_sub.add_argument('--dir','-d',default=".",help='Storage directory')
    parser_sub.add_argument('--txt',action="store_true",help="Add text output")
    parser_sub.add_argument('--csv',action="store_true",help="Add CSV output")
    parser_sub.add_argument('--pdf',action="store_true",help="Add PDF output. This requires pdflatex to be installed")
    parser_sub.add_argument('--add_columns',default=None,type=str,help="Add additional columns found in the mutation database to the text and pdf results")
    parser_sub.add_argument('--meta',default=None,type=str,help="Add meta data from a CSV file to the results. The CSV file must contain a column labelled \"id\" with the same value as the prefix argument")
    parser_sub.add_argument('--verbose','-v',default=0, choices=[0,1,2],help="Verbosity increases from 0 to 2",type=int)
    parser_sub.add_argument('--no_flagstat',action="store_true",help="Don't collect flagstats")
    parser_sub.add_argument('--no_delly',action="store_true",help="Don't collect flagstats")
    parser_sub.add_argument('--delly_bcf_file',help="Precomputed delly bcf file")
    parser_sub.add_argument('--no_coverage',action="store_true",help="Don't collect coverage stats")
    parser_sub.add_argument('--version', action='version', version="TBProfiler version %s" % tbp._VERSION)
    parser_sub.add_argument('--dump_tmp_results',action="store_true",help="Dump temp results")
    parser_sub.add_argument('--coverage_fraction_threshold',default=0,type=int,help='Cutoff used to calculate fraction of region covered by <= this value')
    parser_sub.add_argument('--missing_cov_threshold',default=10,type=int,help='Cutoff used to positions/codons in genes which are missing')
    parser_sub.set_defaults(func=main_profile)


    args = parser.parse_args()
    if vars(args)=={}:
        parser.print_help(sys.stderr)
    else:
        args.func(args)
