#!/usr/bin/env python3
"""
Script for extracting sequence names and lengths from a FASTA file
Argument1: path to FASTA file
Output: csv table printed to STDOUT. First column: sequence names. Second column: sequence lengths.
"""

import argparse
import general_purpose_functions as gpf

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("path", help="Input FASTA file path", type=str)
parser.add_argument("--tab", default=False, action="store_true" , help="Tab separated output instead of comma separated")

args = parser.parse_args()
fasta_data = gpf.read_fasta_in_chunks(args.path)
for header, seq in fasta_data:
    if args.tab == False:
        print(header + "," + str(len(seq)))
    else:
        print(header + "\t" + str(len(seq)))
                
