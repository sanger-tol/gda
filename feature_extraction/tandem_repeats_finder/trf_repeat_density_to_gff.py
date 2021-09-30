#!/usr/bin/env python3
"""
Script for converting Tandem Repeats Finder repeat density values (that have been divided into repeat-rich and repeat-poor regions) into GFF or BED format
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

import sys
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import pandas as pd
import argparse
from genome_decomp_pipeline_shared_functions import print_gff_header_based_on_fasta


def load_repeats_data_as_df(repeat_density_path):
    """
    Input: path to an output file of trf_repeat_density_sliding_window.py
    Output: Input file loaded as a pandas data frame, with a new column segment_classif that gives an ID for each contiguous repeat-rich or repeat-poor region
    """
    previous_rep_rich_flag = None
    previous_scaff = None
    segment_counter = 0
    df = pd.read_csv(repeat_density_path, sep="\t")
    segment_classif = list()
    for row_tuple in df.iterrows():
        row = row_tuple[1]
        rr_flag = row["repeat_rich_flag"]
        scaff = row["scaffold"]
        if rr_flag != previous_rep_rich_flag or scaff != previous_scaff:
            segment_counter += 1
            previous_rep_rich_flag = rr_flag
            previous_scaff = scaff
        segment_classif.append(segment_counter)

    df["segment_classif"] = segment_classif
    return df


def print_reformatted_df(df, out_format, assembly_path):
    """
    Input: 1) output of load_repeats_data_as_df function, 2) output format, either GFF or BED, 3) path to assembly FASTA file, ony needed for GFF output
    Output: either a GFF or a BED file that marks where repeat-rich and repeat-poor regions are in the genome
    """
    region_bed_colours = {"repeat_rich_region": "0,255,0", "repeat_poor_region": "0,0,255"}
    if out_format == "GFF":
        print_gff_header_based_on_fasta(assembly_path)

    segment_count = len(df.segment_classif.unique())
    for segment_i in range(1, segment_count + 1):
        df_chunk = df.loc[df["segment_classif"] == segment_i]
        chunk_start = min(df_chunk["start_pos"]) + 1
        chunk_end = max(df_chunk["end_pos"]) + 1
        chunk_N_fraction = df_chunk["N_fraction"].mean()

        repeat_rich_flag = df_chunk["repeat_rich_flag"].values[0]
        region_type = None
        if repeat_rich_flag == True:
            region_type = "repeat_rich_region"
        elif repeat_rich_flag == False:
            region_type = "repeat_poor_region"

        scaff_name = df_chunk["scaffold"].values[0]
        if out_format == "GFF":
            out_line = scaff_name + "\t" + "TandemRepeatsFinder" + "\t" + region_type + "\t" + str(chunk_start) + "\t" + str(chunk_end) + "\t" + str(chunk_N_fraction) + "\t" + "+" + "\t" + "." + "\t" + "."
        elif out_format == "BED":
            out_line = scaff_name + "\t" + str(chunk_start - 1) + "\t" + str(chunk_end) + "\t" + region_type + "\t" + str(chunk_N_fraction) + "\t" + "+" + "\t" + str(chunk_start - 1) + "\t" + str(chunk_end) + "\t" + region_bed_colours[region_type]
        print(out_line)


def main(repeat_density_path, out_format, assembly_path):
    df = load_repeats_data_as_df(repeat_density_path)
    print_reformatted_df(df, out_format, assembly_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repeat_density_path", type=str, help="Path to an output file of trf_repeat_density_sliding_window.py")
    parser.add_argument("out_format", type=str, help="Output file format (GFF or BED)", choices=["GFF", "BED"])
    parser.add_argument("--assembly_path", type=str, help="Path to genome assembly FASTA file, required if the output is a GFF file")
    args = parser.parse_args()
    if args.out_format == "GFF" and args.assembly_path == None:
        parser.error("GFF option requires the --assembly_path argument")
    main(args.repeat_density_path, args.out_format, args.assembly_path)




