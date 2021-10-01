#!/usr/bin/env python3
"""
Script for running DustMasker to detect low complexity regions in an assembly
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

def fasta_to_uppercase(in_fasta_path, out_fasta_path):
    """
    Converts all sequences in FASTA to uppercase
    """
    fasta_data = gpf.read_fasta_in_chunks(in_fasta_path)
    out_list = list()
    for header, seq in fasta_data:
        out_list.append(">" + header)
        seq = seq.upper()
        seq_list = gpf.split_with_fixed_row_length(seq, 80)
        out_list.extend(seq_list)
    gpf.export_list_as_line_break_separated_file(out_list, out_fasta_path)


def main(in_fasta_path, dustmasker_folder, pipeline_output_folder, chunk_size):
    gpf.run_system_command("mkdir -p " + dustmasker_folder)
    gpf.run_system_command("mkdir -p " + pipeline_output_folder)
    in_fasta_basename = in_fasta_path.split("/")[-1]
    in_fasta_basename = in_fasta_basename.split(".fa")[0]
    uppercase_fasta_path = dustmasker_folder + "/" + in_fasta_basename + "_uppercase.fa"
    fasta_to_uppercase(in_fasta_path, uppercase_fasta_path)
    dustmasker_fasta_path = dustmasker_folder + "/" + in_fasta_basename + "_uppercase_dustmasker.fa"
    bedgraph_path = pipeline_output_folder + "/" + in_fasta_basename + "_dustmasker_low_complexity.bedgraph"

    dustmasker_command = "dustmasker -in {} -out {} -outfmt fasta".format(uppercase_fasta_path, dustmasker_fasta_path)
    gpf.run_system_command(dustmasker_command)
    dustmasker_to_bedgraph_command = "dustmasker_get_masked_seq_percentage_bedgraph.py {} --chunk_size {} > {}".format(dustmasker_fasta_path, chunk_size, bedgraph_path)
    gpf.run_system_command(dustmasker_to_bedgraph_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to a FASTA file where low complexity sequences have been masked using Dustmasker. All nucleotides have to be uppercase in the FASTA file before masking with Dustmasker")
    parser.add_argument("dustmasker_folder", type=str, help="Path to folder for Dustmasker files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph file")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk length (bp) for the estimation of density of low complexity regions using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.in_path, args.dustmasker_folder, args.pipeline_output_folder, args.chunk_size)
