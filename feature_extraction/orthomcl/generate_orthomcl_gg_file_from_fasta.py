#!/usr/bin/env python3
"""
Script for generating gg_file for OrthoMCL from protein FASTA files
Argument1: path to a CSV file. First column: species identifiers (short names). Second column: names of FASTA files for each species (without folder path)
Output: gg_file for OrthoMCL
Argument2: path to folder with protein FASTA files
"""

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


