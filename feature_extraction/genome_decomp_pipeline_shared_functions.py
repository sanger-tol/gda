#!/usr/bin/env python3 
"""
File for functions that are shared between scripts of the decomposition pipeline
"""

import general_purpose_functions as gpf
from collections import OrderedDict

def print_gff_header_based_on_fasta(assembly_path):
    """
    Prints GFF header that contains the scaffold lengths, determined by reading the FASTA file
    """
    print("##gff-version 3")
    fasta_data = gpf.read_fasta_in_chunks(assembly_path)
    for header, seq in fasta_data:
        seq_len = len(seq)
        header = header.split()[0]
        print("##sequence-region" + "   " + header + " 1 " + str(seq_len))


def get_fasta_sequence_lengths(fasta_path):
    """
    Input: path to a FASTA file
    Output: a dictionary where the keys are FASTA headers and values are sequence lengths
    """
    fasta_lengths_dict = OrderedDict()
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    for header, seq in fasta_data:
        header = header.split()[0]
        fasta_lengths_dict[header] = len(seq)
    return fasta_lengths_dict