#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import logging
import os
import sys


# Argument definition function
def parse_arguments():
    parser = argparse.ArgumentParser(description="Filtering transcript assembled by Stringtie", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action='store_true', default=False, help="""Print the  documentation of this tool""")
    parser.add_argument("-bdg", "--bdg_data", default="stdin", help="""Output data from intersect beetwen bedgraph and bed files """)
    parser.add_argument("-o", "--output_gene", default="stdin", help="""output of liste of selected gene""")

    return parser.parse_args()


if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.bdg_data == "stdin":
        inputdata = sys.stdin
    else:
        inputdata = args.bdg_data
        logging.info("Using " + str(inputdata) + "as input data")


    if args.output_gene == "stdin":
        outputdata = sys.stdin
    else:
        outputdata = args.output_gene
        logging.info("Using " + str(outputdata) + "as input data")


    if args.documentation:
        print(str(parse_arguments) + "\n" + str(parse_arguments__doc__))


    # Load Data
    Colnames = [
        'chr_bedgraph', 'start_bedgraph', 'end_bedgraph', 'coverage_bedgraph',
        'chr_bed', 'start_bed', 'end_bed', 'gene_id', 'length', 'strand'
    ]
    bed_data = pd.read_csv(inputdata, sep="\t", names=Colnames)

    # 25th percentile compute by gene_id
    percentiles = (
        bed_data.groupby('gene_id')['coverage_bedgraph']
        .quantile(0.25)
        .rename('p25')
    )
    bed_data = bed_data.merge(percentiles, on='gene_id')

    # Filter line with coverage >= 50th percentile
    filtered_data = bed_data[bed_data['coverage_bedgraph'] >= bed_data['p25']]


    # compute length for each gene_id
    lengths = (
        filtered_data.groupby('gene_id')
        .agg(start_min=('start_bedgraph', 'min'), end_max=('end_bedgraph', 'max'))
        .assign(length=lambda df: df['end_max'] - df['start_min'])
    )

    # Filter gene with length >= 50 and save results
    results = lengths[lengths['length'] >= 50].reset_index()
    results[['gene_id', 'length']].to_csv(outputdata, sep='\t', index=False)

