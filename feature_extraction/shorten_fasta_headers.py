#!/usr/bin/env python3
"""
Script for shortening FASTA headers, by splitting the header and keeping only the first element
"""
# MIT License
# 
# Copyright (c) 2020-2021 Genome Research Ltd.
# 
# Author: Eerik Aunin (ea10@sanger.ac.uk)
# 
# This file is a part of the Genome Decomposition Analysis (GDA) pipeline.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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


