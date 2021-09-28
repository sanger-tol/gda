#!/usr/bin/env python3
"""
Script for shortening FASTA headers, by splitting the header and keeping only the first element
"""

import general_purpose_functions as gpf
import argparse
import sys


def main(fasta_path, delimiter, allow_duplicate_headers):
    headers_list = list()
    in_data = gpf.ll(fasta_path)

    for line in in_data:
        out_line = line
        if line.startswith(">"):
            if delimiter == "":         
                out_line = line.split()[0]
            else:
                out_line = line.split(delimiter)[0]
            if out_line in headers_list and allow_duplicate_headers == False:
                sys.stderr.write("Duplicate FASTA headers ({}) were found in the input file ({}) after truncating the headers with a delimiter\n".format(out_line[1:len(out_line)], fasta_path))
                sys.exit(1)
            headers_list.append(out_line)
        print(out_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("--delimiter", type=str, help="Delimiter string for splitting FASTA headers. Default: any whitespace character", default="")
    parser.add_argument("--allow_duplicate_headers", dest="allow_duplicate_headers", action="store_true")
    args = parser.parse_args()
    main(args.fasta_path, args.delimiter, args.allow_duplicate_headers)

            
