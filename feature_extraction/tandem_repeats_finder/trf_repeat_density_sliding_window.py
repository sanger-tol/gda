#!/usr/bin/env python3
"""
Script for finding repeat density in a genome using sliding window on genome FASTA file where repeats have been masked with Tandem Repeats Finder
Output: tab separated table. Column1: scaffold name. Column2: chunk start coordinate in the scaffold (1-based). Column3: chunk end coordinate in the scaffold.
    Column4: fraction of nucleotides in the chunk that were masked by TandemRepeatsFinder. Column5: True if the fraction of masked nucleotides exceeds a cutoff, False if not.
"""
import sys
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse


def main(in_path, chunk_size, seq_len_cutoff, repeat_rich_cutoff, zero_based):
    fasta_data = gpf.read_fasta_in_chunks(in_path)

    print("scaffold\tstart_pos\tend_pos\tN_fraction\trepeat_rich_flag")
    coords_modifier = int(zero_based)

    for header, seq in fasta_data:
        if len(seq) > seq_len_cutoff:
            seq = seq.upper()
            seq_chunks = gpf.string_to_chunks(seq, chunk_size)
            pos = 1
            for chunk in seq_chunks:
                chunk_N_fraction = chunk.count("N") / len(chunk)
                chunk_repeat_rich_flag = False
                if chunk_N_fraction >= repeat_rich_cutoff:
                    chunk_repeat_rich_flag = True
                chunk_end = pos + len(chunk) - 1
                out_line = "\t".join([header, str(pos - coords_modifier), str(chunk_end - coords_modifier), str(chunk_N_fraction), str(chunk_repeat_rich_flag)])
                print(out_line)
                pos += chunk_size

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to DNA FASTA file where repeats have been masked with Tandem Repeats Finder")
    parser.add_argument("--chunk_size", type=int, help="Sliding window step size (bp)", default=1000)
    parser.add_argument("--seq_len_cutoff", type=int, help="Minimum scaffold size (bp) to be included in the tandem repeat density assessment", default=5000)
    parser.add_argument("--repeat_rich_cutoff", type=float, help="Minimum fraction of tandem repeats a scaffold needs to contain to be considered repeat-rich", default=0.05)
    parser.add_argument("--zero_based", dest="zero_based", action="store_true", help="Use zero-based coordinates instead of 1-based (optional)")
    args = parser.parse_args()
    main(args.in_path, args.chunk_size, args.seq_len_cutoff, args.repeat_rich_cutoff, args.zero_based) 

   