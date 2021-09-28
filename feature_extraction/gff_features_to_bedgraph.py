#!/usr/bin/env python3
"""
Script for converting GFF features to bedgraph format. Output (STDOUT): a bedgraph file with the fractions of regions that contain the query feature in fixed length chunks of scaffolds
"""

import argparse
from genome_decomp_pipeline_shared_functions import get_fasta_sequence_lengths
import general_purpose_functions as gpf
from collections import OrderedDict
import numpy as np


def get_seq_len_dict(gff_data):
    """
    Input: GFF file, loaded as a list
    Output: dictionary where the keys are scaffold names and values are scaffold lengths
    """
    region_data = [n for n in gff_data if n.startswith("##sequence-region")]
    seq_len_dict = OrderedDict()
    for line in region_data:
        split_line = line.split()
        seq_len_dict[split_line[1]] = int(split_line[3])
    return seq_len_dict


def get_blank_feature_arrays(gff_data, assembly_fasta_path):
    """
    Input: GFF file, loaded as a list
    Output: a dictionary where keys are scaffold names and values are 1D boolean numpy arrays (one element for each nucleotide in the scaffold)
    """
    seq_len_dict = None
    if assembly_fasta_path == "":
        seq_len_dict = get_seq_len_dict(gff_data)
    else:
        seq_len_dict = get_fasta_sequence_lengths(assembly_fasta_path)
    arrays_dict = dict()
    
    for selected_scaff in seq_len_dict.keys():
        selected_scaff_len = seq_len_dict[selected_scaff]
        scaff_array = np.zeros(selected_scaff_len, dtype=bool)
        arrays_dict[selected_scaff] = scaff_array
    return arrays_dict


def main(gff_path, query_feature, header_title, assembly_fasta_path, chunk_size):
    gff_data = gpf.l(gff_path)

    arrays_dict = get_blank_feature_arrays(gff_data, assembly_fasta_path)

    for line in gff_data:
        if line.startswith("#") == False:
            split_line = line.split()
            if len(split_line) > 5:
                if split_line[2] == query_feature:
                    feature_scaff = split_line[0]
                    feature_start = int(split_line[3]) - 1
                    feature_end = int(split_line[4])
                    arrays_dict[feature_scaff][feature_start:feature_end] = True

    bedgraph_header = 'track type=bedGraph name="' + header_title + '" description="' + header_title + '" visibility=full color=0,0,255 altColor=0,100,200 priority=20'
    print(bedgraph_header)

    for selected_scaff in arrays_dict:
        scaff_array = arrays_dict[selected_scaff]

        for scaff_start_pos in range(0, scaff_array.size, chunk_size):
            scaff_subset = scaff_array[scaff_start_pos: scaff_start_pos + chunk_size]
            scaff_end_pos = scaff_start_pos + scaff_subset.size
            chunk_mean = np.mean(scaff_subset)
            print(selected_scaff + " " + str(scaff_start_pos) + " " + str(scaff_end_pos) + " " + str(chunk_mean))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff_path", type=str, help="Path to GFF file")
    parser.add_argument("query_feature", type=str, help="Query feature type to search for in the GFF file, e.g. mRNA")
    parser.add_argument("header_title", type=str, help="Title for bedgraph file header")
    parser.add_argument("--assembly_fasta_path", type=str, default="", help="Optional: path to assembly FASTA file. If this is used, assembly sequence lengths are determined by using the FASTA file instead of reading them from the GFF header")
    parser.add_argument("--chunk_size", type=int, help="Chunk size (sliding window step size) in base pairs. Default: 5000", default=5000)
    args = parser.parse_args()
    main(args.gff_path, args.query_feature, args.header_title, args.assembly_fasta_path, args.chunk_size)








