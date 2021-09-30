#!/usr/bin/env python3
"""
File for functions that are shared between scripts of the decomposition pipeline
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
