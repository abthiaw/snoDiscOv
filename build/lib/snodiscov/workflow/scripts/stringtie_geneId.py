#!/usr/bin/env python3

# coding: utf-8

__author__="Thiaw Alphonse"

# Library import
import argparse
import sys
import logging


# Argument definition function
def parse_arguments():
    parser = argparse.ArgumentParser(description="Ce programme permet de modifier les gene_id dans les fichier gtf de stringtie", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-i", "--input", default="stdin", help="""gtf stringtie file""")
    parser.add_argument("-s", "--suffixe", default="stdin", help="""suffixe that will be added to gene_id name""")
    
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.input == "stdin":
        data = sys.stdin
    else:
        data = args.input
        logging.info("Using " + str(data) + "as input data")


    if args.suffixe == "stdin":
        Id_suff = sys.stdin
    else:
        Id_suff = args.suffixe
        logging.info("Using " + str(Id_suff) + "as input data")

    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))
    
    
    with open(data, 'r') as file:
        for li in file:
            li=li.replace("STRG", "STRG"+str(Id_suff))
            print(li, end="")


