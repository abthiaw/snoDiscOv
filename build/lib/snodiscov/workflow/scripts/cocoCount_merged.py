#!/usr/bin/env python3

# conding: utf8
__author__ = "Thiaw Alphonse"

#Library import
import argparse
import logging
import pandas as pd
import os
import sys
from pathlib import Path


#Argument definition function
def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge tpm value of several sample obtained from coco count tools ", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-i", "--pathdir", default="stdin", help="""Path to coco count output directory""")
    parser.add_argument("-f", "--file_out", default="stdin", help="""output file type either tpm  or count file""")
    parser.add_argument("-o", "--tpmvalue", default="stdin", help="""tmp merged dataframe output file name""")

    return parser.parse_args()


if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.pathdir == "stdin":
        inputpath = sys.stdin
    else:
        inputpath = args.pathdir
        logging.info("Using " + str(inputpath) + "as input data")


    if args.tpmvalue == "stdin":
        output = sys.stdin
    else:
        output = args.tpmvalue
        logging.info("Using " + str(output) + "as input data")

    if args.file_out == "stdin":
        File = sys.stdin
    else:
        File = args.file_out
        logging.info("Using " + str(File) + "as input data")


    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))


dcount=pd.DataFrame()
dtpm=pd.DataFrame()

filelist=sorted(Path(inputpath).glob("*.tsv"))  # sorted: deterministic column order
for i, path in enumerate(filelist):
    basename=os.path.basename(path).split('.')[0]
    df=pd.read_csv(path, sep="\t")
    #snodf=df[df.gene_id.str.startswith("STRG")]
    if i==0:
        dfidcount=df[["gene_id", "count"]].rename(columns={"count":basename})
        dfidtpm=df[["gene_id", "tpm"]].rename(columns={"tpm":basename})
    else:
        snodfcount=df[["count"]].rename(columns={"count":basename})
        dcount=pd.concat([dcount,snodfcount], axis=1)
        snodftpm=df[["tpm"]].rename(columns={"tpm":basename})
        dtpm=pd.concat([dtpm,snodftpm], axis=1)

Type=str(File)
if Type=="count":
    dcount=pd.concat([dfidcount,dcount], axis=1)
    dcount=dcount.drop_duplicates()
#print(d)
    dcount.to_csv(output, index=False, sep="\t")
if Type=="tpm":
    dtpm=pd.concat([dfidtpm,dtpm], axis=1)
    dtpm=dtpm.drop_duplicates()
#print(d)
    dtpm.to_csv(output, index=False, sep="\t")

                   
