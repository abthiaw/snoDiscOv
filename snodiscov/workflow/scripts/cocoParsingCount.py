#!/usr/bin/env python3

# coding: utf-8

__author__ = "Thiaw Alphonse Birane"

# Library import
from collections import defaultdict
import argparse
import logging
from pathlib import Path
import os
import sys
import re


# Argument definition function
def parse_arguments():
    parser = argparse.ArgumentParser(description="This program allows to retrieve transcript with 'i' or 'u' and TPM>1 according to CoCo", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-p", "--folder", default="stdin", help="""Path folder of CoCocount output""")
    parser.add_argument("-i", "--gtffile", default="stdin", help="""annotation file of Merged pics""")
    parser.add_argument("-r", "--run_mode", default="stdin", help="""what kind of file do you want to output: annotation file if you make gtfFile as option or transcript coordinate file for bedfile build if make bedFile as option)""")

    return parser.parse_args()


if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.folder == "stdin":
        inputpath = sys.stdin
    else:
        inputpath = args.folder
        logging.info("Using " + str(inputpath) + "as input data")


    if args.gtffile == "stdin":
        inputData = sys.stdin
    else:
        inputData = args.gtffile
        logging.info("Using " + str(inputData) + "as input data")

    if args.run_mode == "stdin":
        mode = sys.stdin
    else:
        mode = args.run_mode
        logging.info("Using " + str(mode) + "as input data")


    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))

    

    dicoCoCo=defaultdict(str)
    repo=inputpath
    files=Path(repo).glob('*.tsv')
    for path in sorted(files):
        with open(path, 'r') as file:
            # snoDiscov v1.0: locate the TPM column from the CoCo header and split
            # on tabs. The former split()[4] collapsed the empty gene_name field
            # of novel transcripts and could read the wrong column (or fail).
            header = file.readline().rstrip('\n').split('\t')
            itpm = header.index('tpm') if 'tpm' in header else 4
            for li in file:
                if li.startswith('MSTRG') or li.startswith('STRG'):
                    li=li.rstrip('\n').split('\t')
                    if float(li[itpm])>=1:
                        dicoCoCo[li[0]]=li[itpm]
                        #print(li[0], li[4], sep='\t')
    

    run=str(mode)
    with open(inputData, 'r') as annotFile:
        for li in annotFile:
            #la=li.rstrip('\t').split()
            la=li.rstrip('\n').split()
            id=re.sub("\"|;","",la[9])
            #print(id)
            if id in dicoCoCo.keys() and run=="gtffile":
                print(li, end="")
            elif id in dicoCoCo.keys() and run=="bedfile":
                if la[2]=="transcript":
                    length=int(la[4])-int(la[3]); strand=la[6]
                    print(la[0],la[3],la[4],id,length,strand,  sep='\t')



