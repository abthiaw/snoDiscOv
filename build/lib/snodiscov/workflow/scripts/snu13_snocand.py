#!/usr/bin/env python3

import argparse
import logging
import pandas as pd
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter expressed transcript (>1TM) in at least 1 sample", formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='-', add_help=True)
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="""My doc string here""")
    parser.add_argument("-doc", "--documentation", action="store_true", default=False, help="""Print the documentation of this tools""")
    parser.add_argument("-i", "--input_count", default="stdin", help="Count file from CoCo count output")
    parser.add_argument("-o", "--output_count", default="stdin", help="output file of expressed (> 1TPM) gene tab ")
    parser.add_argument("-p", "--output_parsing", default="stdin", help="Output for parsing snu13 RIP-PEN-Seq data analysis")

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()

    if args.verbose:
        logging.basicConfig(format='%(asctime)s - %(name)s : %(process)d - %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S', level=logging.INFO)
        logging.info('All information message will get logged')

    if args.input_count == "stdin":
        input_data = sys.stdin
    else:
        input_data = args.input_count
        logging.info("Using" + str(input_data) + "as input data")

    if args.output_count == "stdin":
        output_data = sys.stdin
    else:
        output_data = args.output_count
        logging.info("Using " + str(output_data) + "as input data")    
    
    if args.output_parsing == "stdin":
        output_parsing_data = sys.stdin
    else:
        output_parsing_data = args.output_parsing
        logging.info("Using " + str(output_parsing_data) + "as input data") 

    
    count_df = pd.read_csv(input_data, sep='\t', index_col='gene_id')
    count_snocand = count_df[count_df.index.str.contains('STRG')]
    exp_snocand = count_snocand[(count_snocand > 1).any(axis=1)].reset_index()
    exp_snocand_pasing = exp_snocand[['gene_id']]
    exp_snocand_pasing['SNU13_ripenseq'] = "cd"
    exp_snocand.to_csv(output_data, sep='\t', index=False)
    exp_snocand_pasing.to_csv(output_parsing_data, sep='\t', index=False)
