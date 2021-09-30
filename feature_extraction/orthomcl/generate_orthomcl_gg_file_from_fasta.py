#!/usr/bin/env python3
"""
Script for generating gg_file for OrthoMCL from protein FASTA files
Argument1: path to a CSV file. First column: species identifiers (short names). Second column: names of FASTA files for each species (without folder path)
Output: gg_file for OrthoMCL
Argument2: path to folder with protein FASTA files
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

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf

species_to_fasta_table_path = sys.argv[1]
fasta_folder = sys.argv[2]

species_to_fasta_table = gpf.l(species_to_fasta_table_path)
for line in species_to_fasta_table:
    split_line = line.split(",")
    species = split_line[0]
    fasta_file = split_line[1]
    fasta_full_path = None
    if "/" in fasta_file:
        fasta_full_path = fasta_file
    else:
        fasta_full_path = fasta_folder + "/" + fasta_file
    fasta_headers = list()
    fasta_data = gpf.read_fasta_in_chunks(fasta_full_path)
    for item in fasta_data:
        fasta_headers.append(item[0])
    fasta_headers = [n.split()[0] for n in fasta_headers]
    concat_headers = " ".join(fasta_headers)
    out_line = species + ": " + concat_headers
    print(out_line)


