#!/usr/bin/env python3

#Library import
from collections import defaultdict
import argparse
import logging
import sys
import re


# coding: utf-8

__author__ = "Thiaw Alphonse Birane"


def parse_arguments():
    parser = argparse.ArgumentParser(description="This program craate annotation file based on gtf file from stringtie with target predicted transcripts", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-a", "--gtf", default="stdin", help="""annotation gtf file of stringtie output""")
    parser.add_argument("-o", "--output", default="stdin", help="""gtf annotation output file name""")

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')


    if args.gtf == "stdin":
        gtffile = sys.stdin
    else:
        gtffile = args.gtf
        logging.info("Using " + str(gtffile) + "as input data")

    if args.output == "stdin":
        outputfile = sys.stdin
    else:
        outputfile = args.output
        logging.info("Using " + str(outputfile) + "as input data")


    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))

    """
    dicoStar=defaultdict(str)
    dicoEnd=defaultdict(str)
    with open(bedfile, "r") as bedfile:
        for li in bedfile:
            li=li.rstrip("\n").split()
            chr=li[0]; star=li[1]; end=li[2]; id=chr+':'+star
            dicoStar[id]=star
            dicoEnd[id]=end
    """
    path=(outputfile)
    sys.stdout = open(path, "w")
    set2=set()
    starset=set()
    gene="gene"
    transcript="transcript"
    exon="exon"
    genebiotype='gene_biotype "snoRNA";'
    exonatt='exon_number "1"; exon_version "1";'
    with open(gtffile, "r") as gtfile:
        for k, li in enumerate(gtfile):
            la=li.rstrip("\n").split('\t')
            #chr=la[0]; star=la[3]; end=la[4]; id=chr+':'+star
            if la[2] == "transcript": # Ignore transcript with > 1 exon
                print(la[0],la[1],gene,la[3],la[4],la[5],la[6],la[7],str(la[8]+" "+genebiotype), sep='\t')
                print(la[0],la[1],transcript,la[3],la[4],la[5],la[6],la[7],la[8], sep='\t')
                print(la[0],la[1],exon,la[3],la[4],la[5],la[6],la[7],la[8]+" "+exonatt, sep='\t')
    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs
