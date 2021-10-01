#!/usr/bin/env python3
"""
Script for extracting sequence names and lengths from a FASTA file
Argument1: path to FASTA file
Output: csv table printed to STDOUT. First column: sequence names. Second column: sequence lengths.
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

