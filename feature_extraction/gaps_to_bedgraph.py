#!/usr/bin/env python3
"""
Script for detecting Ns in a FASTA file with sliding window and reporting them in bedgraph format
"""

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
