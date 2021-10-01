#!/usr/bin/env python3
"""
Script for splitting a genome assembly into chunks with fixed size
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
import tempfile
import os


def main(input_fasta_path, output_folder, fasta_chunk_size):
    """
    Splits the input assembly FASTA file into smaller chunk files with fixed size
    """
    input_fasta_path = os.path.abspath(input_fasta_path)
    gpf.run_system_command("mkdir -p " + output_folder)
    output_folder = os.path.abspath(output_folder)

    with tempfile.TemporaryDirectory() as tmpdir:

        assembly_fasta_path = tmpdir + "/assembly.fasta"
        gpf.run_system_command("cp {} {}".format(input_fasta_path, assembly_fasta_path))
        gpf.run_system_command("samtools faidx {}".format(assembly_fasta_path))

        seq_len_file_path = tmpdir + "/seq_len_file.txt"
        gpf.run_system_command("get_fasta_sequence_lengths.py --tab {} > {}".format(assembly_fasta_path, seq_len_file_path))
        windows_bed_file_path = tmpdir + "/fasta_windows.bed"

        gpf.run_system_command("bedtools makewindows -g {} -w {} -i srcwinnum > {}".format(seq_len_file_path, fasta_chunk_size, windows_bed_file_path))
        fasta_windows_path = tmpdir + "/fasta_windows.fa"

        gpf.run_system_command("bedtools getfasta -fi {} -fo {} -bed {}".format(assembly_fasta_path, fasta_windows_path, windows_bed_file_path))
        gpf.run_system_command("split_multifasta_with_fixed_nr_of_sequences.py {} {} 1".format(fasta_windows_path, output_folder))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("output_folder", type=str, help="Path to folder for output files")
    parser.add_argument("--fasta_chunk_size", type=int, help="Length of sequence chunks (bp) into which the input FASTA file will be broken (default: 500000)", default=500000)
    args = parser.parse_args()
    main(args.input_fasta_path, args.output_folder, args.fasta_chunk_size)
