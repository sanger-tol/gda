#!/usr/bin/env python3
"""
Script for dividing a multiFASTA file into smaller FASTA files with a fixed number of sequences
Argument1: path to input FASTA file
Argument2: path for output folder
Argument3: number of sequences per output file
Output: original FASTA file split into individual files
"""

import general_purpose_functions as gpf
import sys
import os

in_path = sys.argv[1]
out_folder = sys.argv[2]
nr_of_seq = int(sys.argv[3])

out_dict = dict()
in_file_basename = os.path.basename(in_path)
in_file_basename = os.path.splitext(in_file_basename)[0]
fasta_data = gpf.read_fasta_in_chunks(in_path)
current_key = -1
for counter, fasta_entry in enumerate(fasta_data):
    if counter % nr_of_seq == 0:
        current_key += 1
        out_dict[current_key] = list()
    out_dict[current_key].append(fasta_entry)

file_name_digits_count = len(str(len(out_dict)))
for out_batch in out_dict:
    out_file_path = out_folder + "/" + str(out_batch).zfill(file_name_digits_count) + "_" + in_file_basename + ".fa"
    batch_data = out_dict[out_batch]
    out_list = list()
    for header, seq in batch_data:
        out_list.append(">" + header)
        out_list.extend(gpf.split_with_fixed_row_length(seq, 80))
    gpf.export_list_as_line_break_separated_file(out_list, out_file_path)




