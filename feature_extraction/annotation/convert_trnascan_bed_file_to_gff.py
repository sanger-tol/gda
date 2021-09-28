#!/usr/bin/env python3
"""
Script for converting tRNAscan output BED file to GFF3 format
"""
import argparse

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import genome_decomp_pipeline_shared_functions
from genome_decomp_pipeline_shared_functions import print_gff_header_based_on_fasta


def main(fasta_path, bed_file_path):
    print_gff_header_based_on_fasta(fasta_path)
    bed_file_data = gpf.l(bed_file_path)

    for line in bed_file_data:
        split_line = line.split()
        scaff = split_line[0]
        feature_start = str(int(split_line[1]) + 1)
        feature_end = split_line[2]
        id = split_line[3]
        score = split_line[4]
        strand = split_line[5]
        out_line = "\t".join((scaff, "tRNAscan", "tRNA", feature_start, feature_end, score, strand, ".", "ID=" + id))
        print(out_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("bed_file_path", type=str, help="Path to tRNAscan output BED file")
    args = parser.parse_args()
    main(args.fasta_path, args.bed_file_path)
