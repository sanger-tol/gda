#!/usr/bin/env python3
"""
Script for making bedgraph tracks that are the sum of all simple repeat tracks or sum of all complex repeat tracks
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
import tempfile


def process_repeats_data_file(merged_repeats_file_path, repeat_type):
    """
    Input: path to a TSV file that consists of merged bedgraph files of a repeat type
    Output (STDOUT): a bedgraph file that contains the sum of the repeat type values for every window
    """
    merged_repeats_data = gpf.ll(merged_repeats_file_path)
    for counter, line in enumerate(merged_repeats_data):
        if counter == 0:
            track_description = "sum_of_" + repeat_type
            out_line = 'track type=bedGraph name="{0}" description="{0}" visibility=full color=0,0,255 altColor=0,100,200 priority=20'.format(track_description)
            print(out_line)
        else:
            split_line = line.split("\t")
            window_end = split_line[2]
            split_line[2] = str(int(window_end) + 1)
            descriptor_columns = " ".join([split_line[4], split_line[1], split_line[2]])

            repeat_values = split_line[5:len(split_line)]
            repeat_values = [float(n) for n in repeat_values]
            repeat_values_sum = sum(repeat_values)
            out_line = descriptor_columns + " " + str(repeat_values_sum)
            print(out_line)


def main(in_folder, assembly_fasta_path, assembly_title, repeat_type, chunk_size):
    with tempfile.TemporaryDirectory() as tmpdirname:
        merged_repeats_file_path = tmpdirname + "/" + repeat_type + "_merged.tsv"
        bedgraph_merging_command = "merge_bedgraph_files {} {} {} {} false > {}".format(assembly_fasta_path, in_folder, chunk_size, assembly_title, merged_repeats_file_path)
        gpf.run_system_command(bedgraph_merging_command, verbose=False)
        process_repeats_data_file(merged_repeats_file_path, repeat_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_folder", type=str, help="Path to folder with bedgraph files of a repeat type")
    parser.add_argument("assembly_fasta_path", type=str, help="Path to the FASTA file of the assembly that corresponds to the bedgraph files")
    parser.add_argument("assembly_title", type=str, help="Species name for the assembly")
    parser.add_argument("repeat_type", type=str, help="Repeat type in the input bedgraph files", choices=["simple_repeats", "complex_repeats"])
    parser.add_argument("--chunk_size", type=int, help="Genome chunk length (bp) of sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.in_folder, args.assembly_fasta_path, args.assembly_title, args.repeat_type, args.chunk_size)


