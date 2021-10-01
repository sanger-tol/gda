#!/usr/bin/env python3
"""
Script for detection of ectopic organellar sequences using BLAST
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
import argparse
import os
import sys

def get_organellar_seq_header(ref_organellar_fasta_path):
    """
    Extracts the sequence name from organellar sequence FASTA file header
    """
    ref_organellar_fasta_data = gpf.l(ref_organellar_fasta_path)
    ref_organellar_fasta_header = ref_organellar_fasta_data[0]
    ref_organellar_fasta_header = ref_organellar_fasta_header[1:len(ref_organellar_fasta_header)]
    ref_organellar_fasta_header = ref_organellar_fasta_header.split()[0]
    return ref_organellar_fasta_header


def check_for_ectopic_organellar_seq(organellar_ref_fasta_path, assembly_fasta_path, out_folder, pipeline_output_folder, seq_tag, chunk_size, threads, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff):
    """
    Runs BLAST to try to detect ectopic organellar sequences (either mitochondrion of apicoplast) in the nuclear genome assembly
    """
    blastdb_path = out_folder + "/ref_" + seq_tag + ".fa"
    makeblastdb_stdout_path = out_folder + "/ref_" + seq_tag + "_makeblastdb_stdout.txt"
    assembly_fasta_filename = assembly_fasta_path.split("/")[-1]
    assembly_fasta_basename = assembly_fasta_filename.split(".")[0]

    blast_output_file_path = out_folder + "/" + assembly_fasta_basename + "_" + seq_tag + "_blast.txt"
    ectopic_organellar_gff_path = out_folder + "/" + assembly_fasta_basename + "_ectopic_" + seq_tag + ".gff"
    ectopic_organellar_bedgraph_path = pipeline_output_folder + "/" + assembly_fasta_basename + "_ectopic_" + seq_tag + ".bedgraph"

    ref_organellar_fasta_header = get_organellar_seq_header(organellar_ref_fasta_path)
    gpf.run_system_command("cp {} {}".format(organellar_ref_fasta_path, blastdb_path))
    gpf.run_system_command("makeblastdb -in {} -out {} -dbtype nucl > {}".format(blastdb_path, blastdb_path, makeblastdb_stdout_path))

    blast_command = "blastn -db {} -query {} -out {} -evalue 1e-30 -num_threads {} -outfmt 6".format(blastdb_path, assembly_fasta_path, blast_output_file_path, str(threads))
    gpf.run_system_command(blast_command)
    if os.stat(blast_output_file_path).st_size == 0:
        sys.stderr.write("No BLAST matches to the organellar sequence (" + organellar_ref_fasta_path + ") were found in the assembly\n")
    else:
        gpf.run_system_command("filtered_ectopic_organellar_blast_hits_to_gff.py {} {} {} {} {} --ectopic_organellar_length_cutoff {} --ectopic_organellar_pident_cutoff {}".format(blast_output_file_path, assembly_fasta_path, ectopic_organellar_gff_path, seq_tag, ref_organellar_fasta_header, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff))
        gpf.run_system_command("gff_features_to_bedgraph.py {} {} {} --chunk_size {} > {}".format(ectopic_organellar_gff_path, "putative_ectopic_" + seq_tag, "putative_ectopic_" + seq_tag, chunk_size, ectopic_organellar_bedgraph_path))


def main(assembly_fasta_path, ref_mitoch_fasta_path, out_folder, pipeline_output_folder, ref_apicoplast_fasta_path, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff, chunk_size, threads):
    gpf.run_system_command("mkdir -p " + out_folder)
    gpf.run_system_command("mkdir -p " + pipeline_output_folder)
    if ref_mitoch_fasta_path != "NA":
        check_for_ectopic_organellar_seq(ref_mitoch_fasta_path, assembly_fasta_path, out_folder, pipeline_output_folder, "mitochondrion", chunk_size, threads, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff)
    if ref_apicoplast_fasta_path != "NA":
        check_for_ectopic_organellar_seq(ref_apicoplast_fasta_path, assembly_fasta_path, out_folder, pipeline_output_folder, "apicoplast", chunk_size, threads, ectopic_organellar_length_cutoff, ectopic_organellar_pident_cutoff)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("ref_mitoch_fasta_path", type=str, help="Path to reference mitochondrion FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder path for output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder path for main output files of the pipeline (where bedgraph files will be saved)")
    parser.add_argument("--ref_apicoplast_fasta_path", type=str, help="Optional: path to reference apicoplast FASTA file", default="NA")
    parser.add_argument("--ectopic_organellar_length_cutoff", type=int, help="Minimum BLAST alignment length cutoff for identifying ectopic organellar sequences (default: 150)", default=150)
    parser.add_argument("--ectopic_organellar_pident_cutoff", type=int, help="Minimum BLAST alignment percentage identity cutoff for identifying ectopic organellar sequences (default: 80)", default=80)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--threads", type=int, help="Number of threads for running BLAST (default: 1)", default=1)
    args = parser.parse_args()
    main(args.assembly_fasta_path, args.ref_mitoch_fasta_path, args.out_folder, args.pipeline_output_folder, args.ref_apicoplast_fasta_path, args.ectopic_organellar_length_cutoff, args.ectopic_organellar_pident_cutoff, args.chunk_size, str(args.threads))






