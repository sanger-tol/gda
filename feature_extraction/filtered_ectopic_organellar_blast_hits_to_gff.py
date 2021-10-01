#!/usr/bin/env python3
"""
Script for filtering the BLAST results to detect ectopic mitochondrial and apicoplast sequences. The script outputs the ectopic mitochondrion and apicoplast BLAST hits as GFF
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

import argparse
import sys
import pandas as pd
from collections import defaultdict
from genome_decomp_pipeline_shared_functions import get_fasta_sequence_lengths


def get_putative_organellar_scaffolds(blast_df, fasta_lengths_dict):
    """
    Function for detecting scaffolds that are probably the real mitochondrion or apicoplast (as opposed to ectopic mitochondrion or apicoplast).
        This is determined by the percentage of the scaffold that gets aligned with a reference organelle sequence.
    Input: 1) data frame with BLAST results, 2) dictionary with FASTA sequence lengths
    Output: list of putative organellar scaffolds
    """
    alignment_lengths_dict = defaultdict(int)
    for row_tuple in blast_df.iterrows():
        row = row_tuple[1]
        alignment_lengths_dict[row["qseqid"]] += row["length"]

    putative_organellar_scaffs = list()

    for fasta_header in fasta_lengths_dict:
        seq_len = fasta_lengths_dict[fasta_header]
        if fasta_header in alignment_lengths_dict:
            aligned_seq_len_percentage = (alignment_lengths_dict[fasta_header] / seq_len) * 100
            if aligned_seq_len_percentage >= 90:
                putative_organellar_scaffs.append(str(fasta_header))
    return putative_organellar_scaffs


def filter_blast_hits(df, seq_type, target_seqid, fasta_lengths_dict, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff):
    """
    Filters BLAST hits data frame to remove sequences that are not ectopic or do not pass the ectopic organellar sequence length or percentage identity cutoffs
    """
    df2 = df.loc[df["sseqid"] == target_seqid]

    putative_organellar_scaffs = get_putative_organellar_scaffolds(df2, fasta_lengths_dict)
    if len(putative_organellar_scaffs) > 0:
        sys.stderr.write("Putative non-ectopic organellar scaffold(s) detected (" + seq_type + "): " + ", ".join(putative_organellar_scaffs) + "\n")
    else:
        sys.stderr.write("No putative organellar scaffold(s) were detected\n")

    for putative_organellar_scaff in putative_organellar_scaffs:
        df2 = df2.drop(df2[df2.qseqid == putative_organellar_scaff].index)

    df2 = df2.drop(df2[df2.length < ectopic_organellar_length_cutoff].index)
    df2 = df2.drop(df2[df2.pident < ectopic_organellar_pident_cutoff].index)
    return df2


def output_selected_blast_hits_as_gff(gff_path, df2, fasta_lengths_dict, seq_type, target_seqid):
    """
    Writes the coordinates of the putative ectopic organellar regions as a GFF file
    """
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
        for scaff_name, scaff_len in fasta_lengths_dict.items():
            scaff_name = scaff_name.split()[0]
            f.write("##sequence-region   {} 1 {}\n".format(str(scaff_name), str(scaff_len)))

        for counter, row_tuple in enumerate(df2.iterrows()):
            row = row_tuple[1]
            out_line = "\t".join([row["qseqid"], "blast", "putative_ectopic_" + seq_type, str(row["qstart"]), str(row["qend"]), str(row["bitscore"]), ".", ".", target_seqid + "_match_" + str(counter).zfill(3), "\n"])
            f.write(out_line)


def main(blast_results_path, assembly_fasta_path, gff_path, seq_type, target_seqid, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff):
    fasta_lengths_dict = get_fasta_sequence_lengths(assembly_fasta_path)

    df = pd.read_csv(blast_results_path, sep="\t", header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    df2 = filter_blast_hits(df, seq_type, target_seqid, fasta_lengths_dict, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff)
    output_selected_blast_hits_as_gff(gff_path, df2, fasta_lengths_dict, seq_type, target_seqid)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blast_results_path", type=str, help="Path to BLAST results file (reference organellar sequence vs the assembly, outfmt 6)")
    parser.add_argument("assembly_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("gff_path", type=str, help="Path for output GFF file for ectopic organellar sequences")
    parser.add_argument("seq_type", type=str, help="Organellar sequence type ('mitochondrion' or 'apicoplast')")
    parser.add_argument("target_seqid", type=str, help="Sequence ID (from FASTA header) of the reference organellar sequence")
    parser.add_argument("--ectopic_organellar_length_cutoff", type=int, help="Minimum BLAST alignment length cutoff for identifying ectopic organellar sequences (default: 150)", default=150)
    parser.add_argument("--ectopic_organellar_pident_cutoff", type=int, help="Minimum BLAST alignment percentage identity cutoff for identifying ectopic organellar sequences (default: 80)", default=80)
    args = parser.parse_args()
    main(args.blast_results_path, args.assembly_fasta_path, args.gff_path, args.seq_type, args.target_seqid, args.ectopic_organellar_length_cutoff, args.ectopic_organellar_pident_cutoff)
