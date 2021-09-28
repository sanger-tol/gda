#!/usr/bin/env python3
"""
Script for shortening FASTA headers, by splitting the header and keeping only the first element
"""

import general_purpose_functions as gpf
import argparse

def main(fasta_path, delimiter):
    in_data = gpf.ll(fasta_path)

    for line in in_data:
        out_line = line
        if line.startswith(">"):         
            out_line = line.split(delimiter)[0]
        print(out_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("--delimiter", type=str, help="Delimiter string for splitting FASTA headers. Default: ' ' (space)", default=" ")
    args = parser.parse_args()
    main(args.fasta_path, args.delimiter)

            
