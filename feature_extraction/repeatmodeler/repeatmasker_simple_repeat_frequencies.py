#!/usr/bin/env python3
"""
Script for finding simple repeat frequencies in in GFF derived from the output of RepeatModeler + RepeatMasker
"""

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import pandas as pd
import argparse

def main(in_path, out_path):
    motifs_dict = dict()
    in_data = gpf.l(in_path)
    for line in in_data:
        if line.startswith("#") == False:
            seq = gpf.spl(line, "Motif:(", ")n")
            rev_comp_seq = gpf.reverse_complement(seq)
            if rev_comp_seq in motifs_dict:
                seq = rev_comp_seq

            if seq not in motifs_dict:
                motifs_dict[seq] = {"repeat_units_count": 0, "occurrences": 0, "total_nucleotides": 0, "gc": None}

            split_line = line.split()
            repeat_start = int(split_line[3])
            repeat_end = int(split_line[4])
            seq_len = len(seq)
            repeat_length = repeat_end - repeat_start
            repeat_units_count = repeat_length / seq_len
            motifs_dict[seq]["repeat_units_count"] += repeat_units_count
            motifs_dict[seq]["occurrences"] += 1
            motifs_dict[seq]["total_nucleotides"] += repeat_length
            gc = ((seq.count("G") + seq.count("C")) / seq_len) * 100
            motifs_dict[seq]["gc"] = gc


    df = pd.DataFrame(motifs_dict)
    df = df.transpose()
    df.to_csv(out_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to input GFF file (simple repeats output of condense_simple_repeat_sequences.py)")
    parser.add_argument("out_path", type=str, help="Path for output CSV file")
    args = parser.parse_args()
    main(args.in_path, args.out_path)




