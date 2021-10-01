#!/usr/bin/env python3
"""
Script for converting RepeatMasker repeat coordinates from GFF to bedgraph
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
import argparse
from collections import defaultdict

def get_selected_repeat_families_list(repeats_gff_path, occurrences_cutoff):
    """
    Input: 1) path to repeats GFF from RepeatMasker, reformatted using reformat_repeatmasker_gff.py, 2) minimum number of occurrences for the repeat type to be included in bedgraph output
    Output: names of repeat familes that will be included in bedgraph output
    """
    gff_data = gpf.l(repeats_gff_path)
    gff_data = [n for n in gff_data if n.startswith("#") == False]

    repeat_families_dict = defaultdict(int)
    for line in gff_data:
        split_line = line.split("\t")
        repeat_family_name = split_line[2]
        repeat_families_dict[repeat_family_name] += 1
    selected_repeat_families = list()

    for repeat_family_name, count in repeat_families_dict.items():
        if count >= occurrences_cutoff:
            selected_repeat_families.append(repeat_family_name)
    return selected_repeat_families


def main(repeats_gff_path, out_folder, assembly_prefix, repeat_type_prefix, occurrences_cutoff, chunk_size, assembly_fasta_path):

    gpf.run_system_command("mkdir -p " + out_folder)
    selected_repeat_families = get_selected_repeat_families_list(repeats_gff_path, occurrences_cutoff)

    for selected_repeat_family in selected_repeat_families:
        selected_repeat_family_title = selected_repeat_family
        if selected_repeat_family_title.startswith("("):
            selected_repeat_family_title = selected_repeat_family_title[1:len(selected_repeat_family_title) - 2]
        out_path = out_folder + "/" + assembly_prefix + "_" + repeat_type_prefix + "_" + selected_repeat_family_title + ".bedgraph"
        bedgraph_track_title = repeat_type_prefix + "_" + selected_repeat_family_title
        os_command = "gff_features_to_bedgraph.py \'{}\' \'{}\' \'{}\' --chunk_size {} --assembly_fasta_path \'{}\' > {}".format(repeats_gff_path, selected_repeat_family, bedgraph_track_title, str(chunk_size), assembly_fasta_path, out_path)
        gpf.run_system_command(os_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repeats_gff_path", type=str, help="path to repeats GFF from RepeatMasker, reformatted using reformat_repeatmasker_gff.py")
    parser.add_argument("out_folder", type=str, help="Path to folder for output files")
    parser.add_argument("assembly_prefix", type=str, help="Identifier for the assembly (for use in output file names)")
    parser.add_argument("repeat_type_prefix", type=str, help="Repeat type identifier", choices=["simple_repeats", "complex_repeats"])
    parser.add_argument("--occurrences_cutoff", type=int, help="Minimum number of occurrences for a repeat type to be included in bedgraph output (default: 20)", default=20)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--assembly_fasta_path", type=str, default="", help="Optional: path to assembly FASTA file. If this is used, assembly sequence lengths are determined by using the FASTA file instead of reading them from the GFF header")
    args = parser.parse_args()
    main(args.repeats_gff_path, args.out_folder, args.assembly_prefix, args.repeat_type_prefix, args.occurrences_cutoff, args.chunk_size, args.assembly_fasta_path)


