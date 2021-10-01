#!/usr/bin/env python3
"""
Script for converting the output of EMBOSS einverted to GFF3 format
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
from genome_decomp_pipeline_shared_functions import print_gff_header_based_on_fasta


def print_coords_collection_as_gff(coords_collection, fasta_path):
    """
    Prints the GFF3 output, which includes the scaffold lengths from the FASTA file and the coordinates of inverted repeats
    """
    print_gff_header_based_on_fasta(fasta_path)
    ids_digits_count = len(str(len(coords_collection)))
    for counter, repeat_dict in enumerate(coords_collection):
        repeat_id = str(counter).zfill(ids_digits_count)
        out_list = [repeat_dict["scaff_name"], "einverted", "inverted_repeat", repeat_dict["start_coord1"], repeat_dict["end_coord2"], repeat_dict["score"], "+", "0", "ID=inverted_repeat_" + repeat_id]
        print("\t".join(out_list))
        out_list = [repeat_dict["scaff_name"], "einverted", "inverted_repeat_half", repeat_dict["start_coord1"], repeat_dict["end_coord1"], repeat_dict["score"], "+", "0", "ID=inverted_repeat_" + repeat_id + "_left;Parent=inverted_repeat_" + repeat_id]
        print("\t".join(out_list))
        out_list = [repeat_dict["scaff_name"], "einverted", "inverted_repeat_half", repeat_dict["start_coord2"], repeat_dict["end_coord2"], repeat_dict["score"], "+", "0", "ID=inverted_repeat_" + repeat_id + "_right;Parent=inverted_repeat_" + repeat_id]
        print("\t".join(out_list))


def main(einverted_file_path, fasta_path):
    in_data = gpf.ll(einverted_file_path)
    repeat_dict = None

    coords_collection = list()
    line_counter = 0
    for line in in_data:
        mod_val = line_counter % 5

        if mod_val == 0:
            if repeat_dict is not None:
                coords_collection.append(repeat_dict)
            repeat_dict = {"scaff_name": None, "score": None, "match_percentage": None, "gaps": None, "start_coord1": None, "end_coord1": None, "start_coord2": None, "end_coord2": None}
        elif mod_val == 1:
            repeat_dict["scaff_name"] = line.split(":")[0]
            line = ": Score" + line.split(": Score")[1]
            repeat_dict["score"] = gpf.spl(line, "Score ", ":")
            match_percentage = gpf.spl(line, "(", "%)")
            match_percentage = match_percentage.strip()
            repeat_dict["match_percentage"] = match_percentage
            gaps = gpf.spl(line, ", ", "gap")
            gaps = gaps.strip()
            repeat_dict["gaps"] = gaps
        elif mod_val == 2:
            split_line = line.split()
            repeat_dict["start_coord1"] = split_line[0]
            repeat_dict["end_coord1"] = split_line[-1]
        elif mod_val == 4:
            split_line = line.split()
            repeat_dict["start_coord2"] = split_line[-1]
            repeat_dict["end_coord2"] = split_line[0]

        line_counter += 1

    if repeat_dict is not None:
        coords_collection.append(repeat_dict)
    print_coords_collection_as_gff(coords_collection, fasta_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("einverted_file_path", type=str, help="Path to input *.einverted file (output of EMBOSS einverted)")
    parser.add_argument("fasta_path", type=str, help="Path to FASTA file that was used to generate the *.einverted file")
    args = parser.parse_args()
    main(args.einverted_file_path, args.fasta_path)


