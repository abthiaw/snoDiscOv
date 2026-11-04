#!/usr/bin/env python3

""" Your Name Code (version 1.0.0) """

__author__ = "Thiaw Alphonse Birane"

#Library import
from collections import defaultdict
import argparse
import logging
import re
import sys

#Argument definition function
def parse_arguments():
    parser = argparse.ArgumentParser(description=" This program goal to generate annotation file or transcripts coordonates tab file predicted by stringtie and have class code *i* or *u* according to gffcompare tool, transcripts size between 50 and 200nt was kept", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-i", "--inputFile1", default="stdin", help=""".tmap file of gffcompare output""")
    parser.add_argument("-j", "--inputFile2", default="stdin", help="""annotation gtf file of stringtie output""")
    parser.add_argument("-r", "--run_option", default="stdin", help="""what kind of file do you want to output: annotation file if you make annotation as option or transcript coordinate file for bedfile build if make coordinate as option)""")

    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.inputFile1 == "stdin":
        inputData1 = sys.stdin
    else:
        inputData1 = args.inputFile1
        logging.info("Using " + str(inputData1) + "as input data")
    

    if args.inputFile2 == "stdin":
        inputData2 = sys.stdin
    else:
        inputData2 = args.inputFile2
        logging.info("Using " + str(inputData2) + "as input data")

    if args.run_option == "stdin":
        runing_option = sys.stdin
    else:
        runing_option = args.run_option
        logging.info("Using " + str(runing_option) + "as input data")

    
    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))   



    dico=defaultdict(str)
    with open(inputData1, "r") as gffcf:
        title = gffcf.readline()
        for li in gffcf:
            li=li.rstrip("\n").split("\t")
            #print(li[7])
            if li[2]=="u" or li[2]=="i":
                dico[li[4]]=li[2]


    Torun=str(runing_option)
    with open(inputData2, "r") as annotf:
        for li in (annotf): 
            if not li.startswith("#"):
                la=li.rstrip("\n").split()
                #id=la[11]
                best_id=re.sub(r'"|;','',la[11])
                #tpm=la[12]
                #print(best_id)
                if best_id in dico:
                    length=int(la[4])-int(la[3])
                    if length>=50 and length<=200 and Torun=="coordinate":
                        print(la[0], la[3], la[4], la[6], length, dico[best_id], sep='\t')
                    if length>=50 and length<=200 and Torun=="annotation":
                        print(li, end='')





