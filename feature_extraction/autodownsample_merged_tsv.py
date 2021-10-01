#!/usr/bin/env python3
"""
Script for automatically triggering the downsampling of the TSV table that has been derived by merging bedgraph files
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
import sys
from genome_decomp_pipeline_shared_functions import get_assembly_size


def main(fasta_path, merged_tsv_path, chunk_size, target_number_of_windows):
    if chunk_size <= 0:
        sys.stderr.write("chunk_size cannot be 0 or a negative number\n")
        sys.exit(1)
    if target_number_of_windows <= 0:
        sys.stderr.write("target_number_of_windows cannot be 0 or a negative number\n")
        sys.exit(1)

    assembly_size = get_assembly_size(fasta_path)

    optimal_window_size = assembly_size / target_number_of_windows
    downsampling_factor_float = optimal_window_size / chunk_size
    downsampling_factor = int(round(downsampling_factor_float))

    if downsampling_factor > 1:
        output_file_path = merged_tsv_path.rstrip(".tsv") + "_downsampled_{}x.tsv".format(downsampling_factor)
        downsampling_command = "downsample_merged_bedgraph_file {} {} > {}".format(merged_tsv_path, downsampling_factor, output_file_path)
        gpf.run_system_command(downsampling_command)
    else:
        sys.stderr.write("The TSV table that has been produced does not appear to need downsampling to have {} or fewer windows\n".format(target_number_of_windows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("merged_tsv_path", type=str, help="Path to the TSV file that has been derived by merging bedgraph files")
    parser.add_argument("--chunk_size", type=int, help="Sliding window step size (bp) for running GDA (default: 5000)", default=5000)
    parser.add_argument("--target_number_of_windows", type=int, help="Target number of windows that should be in the TSV after downsampling (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.fasta_path, args.merged_tsv_path, args.chunk_size, args.target_number_of_windows)





