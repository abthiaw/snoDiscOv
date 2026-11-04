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
    parser = argparse.ArgumentParser(description="Ce programme permet de creer un fichier d'annotation a partir du fichier gtf de stringtie et de transcrits cible dans ce fichier sous format bed ", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-b", "--bed", default="stdin", help="""bed file of target transcript in stringtie gtf file""")
    parser.add_argument("-a", "--gtf", default="stdin", help="""annotation gtf file of stringtie output""")
    parser.add_argument("-o", "--output", default="stdin", help="""gtf annotation output file name""")

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.bed == "stdin":
        bedfile = sys.stdin
    else:
        bedfile = args.bed
        logging.info("Using " + str(bedfile) + "as input data")

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


    dicoStar=defaultdict(str)
    dicoEnd=defaultdict(str)
    with open(bedfile, "r") as bedfile:
        for li in bedfile:
            li=li.rstrip("\n").split()
            chr=li[0]; star=li[1]; end=li[2]; id=chr+':'+star
            dicoStar[id]=star
            dicoEnd[id]=end

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
            chr=la[0]; star=la[3]; end=la[4]; id=chr+':'+star
            # snoDiscov v1.0 fix: the first line used to be printed even when its
            # start is not a merged-locus start (id not in dicoStar), producing a
            # gene with empty coordinates. It now goes through the same test as
            # the other lines; output is unchanged whenever the old code was valid.
            if k == 0:
                if id in dicoStar.keys():
                    print(la[0],la[1],gene,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],str(la[8]+" "+genebiotype), sep='\t')
                    print(la[0],la[1],transcript,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8], sep='\t')
                    print(la[0],la[1],exon,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8]+" "+exonatt, sep='\t')
                    starset.add(id)
            else:
                if len(starset) != 0:                      
                    if id in dicoStar.keys() and id not in starset:
                        if la[2]=='transcript':
                            print(la[0],la[1],gene,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],str(la[8]+" "+genebiotype), sep='\t')
                            print(la[0],la[1],transcript,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8], sep='\t')
                            print(la[0],la[1],exon,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8]+" "+exonatt, sep='\t')
                        elif la[2]=='exon':
                            print(la[0],la[1],gene,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],str(la[8]+" "+genebiotype), sep='\t')
                            print(la[0],la[1],transcript,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8], sep='\t')
                            print(la[0],la[1],exon,dicoStar[id],dicoEnd[id],la[5],la[6],la[7],la[8]+" "+exonatt, sep='\t')
                        starset.add(id)
                else:
                    starset.add(id)
    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs
