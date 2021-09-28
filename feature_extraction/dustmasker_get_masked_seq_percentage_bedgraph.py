#!/usr/bin/env python3
"""
Script for finding the percentages of nucleotides that got masked by Dustmasker in a FASTA file, as a sliding window
"""

import general_purpose_functions as gpf
import argparse

def get_masked_seq_percentage(seq):
    """
    Input: sequence where some regions have been masked with Dustmasker
    Output: masked sequence percentage
    """
    lowercase_count = seq.count("a") + seq.count("t") + seq.count("g") + seq.count("c")
    uppercase_count = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
    masked_percentage = "NA"
    total_count = lowercase_count + uppercase_count
    if total_count > 0:
        masked_percentage = (lowercase_count / total_count) * 100
    return masked_percentage


def main(in_path, chunk_size):
    header_title = "dustmasker_low_complexity_percentage"

    bedgraph_header = 'track type=bedGraph name="' + header_title + '" description="' + header_title + '" visibility=full color=0,0,255 altColor=0,100,200 priority=20'
    print(bedgraph_header)

    in_data = gpf.read_fasta_in_chunks(in_path)

    for header, seq in in_data:
        pos = 0
        seq_chunks = gpf.string_to_chunks(seq, chunk_size)
        for seq_chunk in seq_chunks:
            masked_percentage = get_masked_seq_percentage(seq_chunk)
            end_pos = pos + len(seq_chunk)
            if masked_percentage != "NA":
                print(header, pos, end_pos, masked_percentage)
            pos = end_pos


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to a FASTA file where low complexity sequences have been masked using Dustmasker. All nucleotides have to be uppercase in the FASTA file before masking with Dustmasker")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk length (bp) for the estimation of density of low complexity regions using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.in_path, args.chunk_size)

