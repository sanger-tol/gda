#!/usr/bin/env python3
"""
Script for condensing a list of simple repeat sequences to remove redundant sequences.
For example: TAA and TTA are the same sequence, one is the reverse complement of the other. TAATAA repeat is the same as TAA repeat. AATAAT is the same as TAATAA but with a shifted starting point
Input: GFF with repeat locations from RepeatMasker, processed with process_repeatmasker_gffs.py to extract only simple repeat sequences
Output: simple repeats GFF with redundant sequences collapsed into one sequence
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

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse
import re


def get_rotated_string_versions(str):
    #https://www.geeksforgeeks.org/generate-rotations-given-string/
    # Print all the rotated strings.
    # This code is contributed
    # by SHUBHAMSINGH10

    lenn = len(str)

    # Generate all rotations
    # one by one and print
    temp = [0] * (lenn)
    for i in range(lenn):
        j = i # Current index in str
        k = 0 # Current index in temp

        # Copying the second part from
        # the point of rotation.
        while (j < len(str)):

            temp[k] = str[j]
            k += 1
            j += 1


        # Copying the first part from
        # the point of rotation.
        j = 0
        while (j < i) :

            temp[k] = str[j]
            j += 1
            k += 1
        yield "".join(temp)


def repeated(s, REPEATER):
    # https://stackoverflow.com/questions/29481088/how-can-i-tell-if-a-string-repeats-itself-in-python
    match = REPEATER.match(s)
    return match.group(1) if match else None


def gff_line_extract_motif(gff_line):
    """
    Extracts the sequence motif from a line of the input GFF file
    """
    motif = gff_line.split()[2]
    motif = motif[1:len(motif) - 2]
    return motif



def get_motifs_dict(gff_data):
    """
    Function for detecting redundant repeat motifs by checking the reverse complement of the motifs and also all possible rotations of the circular string
    Input: GFF file with simple repeats, loaded as a list of strings
    Output: dictionary. Keys: repeat motifs in the form as they appear in the GFF
        Values: if the key (or its reverse complement or rotation) matches a previously encountered motif, the value will be this previously encountered motif.
            If the key and its reverse complement or rotations do not match any previously encountered motif, the value will be the same as the key.
    """
    gff_data = [n for n in gff_data if n.startswith("#") == False]
    REPEATER = re.compile(r"(.+?)\1+$")
    used_motifs = list()
    motifs_dict = dict()
    for line in gff_data:
        motif = gff_line_extract_motif(line)
        if motif not in motifs_dict:
            collapsed_motif = motif
            sub = repeated(motif, REPEATER)
            if sub:
                collapsed_motif = sub
            rev_comp_collapsed_motif = gpf.reverse_complement(collapsed_motif)
            if collapsed_motif in used_motifs:
                motifs_dict[motif] = collapsed_motif
            elif rev_comp_collapsed_motif in used_motifs:
                motifs_dict[motif] = rev_comp_collapsed_motif
            else:
                rotated_collapsed_motif = list(get_rotated_string_versions(collapsed_motif))
                rotated_collapsed_motif.extend(list(get_rotated_string_versions(rev_comp_collapsed_motif)))
                rotated_collapsed_motif = [n for n in rotated_collapsed_motif if n in used_motifs]
                if len(rotated_collapsed_motif) > 0:
                    motifs_dict[motif] = rotated_collapsed_motif[0]
                else:
                    motifs_dict[motif] = collapsed_motif
                    used_motifs.append(collapsed_motif)
    return motifs_dict


def main(gff_path):
    gff_data = gpf.l(gff_path)
    motifs_dict = get_motifs_dict(gff_data)
    for line in gff_data:
        if line.startswith("#"):
            print(line)
        else:
            motif = gff_line_extract_motif(line)
            output_motif = motifs_dict[motif]
            if output_motif != motif:
                motif = "(" + motif + ")n"
                output_motif = "(" + output_motif + ")n"
                split_line = line.split("\t")
                split_line[2] = output_motif
                out_line = "\t".join(split_line)
                print(out_line)
            else:
                print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff_path", type=str, help="Path to input GFF file")
    args = parser.parse_args()
    main(args.gff_path)
