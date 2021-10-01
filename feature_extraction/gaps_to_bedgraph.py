#!/usr/bin/env python3
"""
Script for detecting Ns in a FASTA file with sliding window and reporting them in bedgraph format
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


def main(fasta_path, chunk_size):
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    print('track type=bedGraph name="N_percentage" description="N_percentage" visibility=full color=0,0,255 altColor=0,100,200 priority=20"')

    for header, seq in fasta_data:
        header = header.split()[0]
        seq = seq.upper()
        split_seq = gpf.string_to_chunks(seq, chunk_size)
        for counter, seq_chunk in enumerate(split_seq):
            chunk_N_count = seq_chunk.count("N")
            chunk_N_percentage = (chunk_N_count / len(seq_chunk)) * 100
            window_coord = counter * chunk_size
            print(header, window_coord, window_coord + len(seq_chunk), chunk_N_percentage)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.fasta_path, args.chunk_size)
