#!/usr/bin/env python3
"""
Script for running einverted to detect inverted repeats
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

def main(assembly_fasta_path, out_folder, pipeline_output_folder, chunk_size):

    assembly_fasta_filename = assembly_fasta_path.split("/")[-1]
    assembly_fasta_basename = assembly_fasta_filename.split(".")[0]

    einverted_output_fasta_path = out_folder + "/" + assembly_fasta_basename + "_einverted.fa"
    einverted_alignments_path = out_folder + "/" + assembly_fasta_basename + ".einverted"
    einverted_gff_path = out_folder + "/" + assembly_fasta_basename + "_einverted.gff"
    einverted_bedgraph_path = pipeline_output_folder + "/" + assembly_fasta_basename + "_einverted.bedgraph"

    gpf.run_system_command("mkdir -p " + out_folder)
    gpf.run_system_command("mkdir -p " + pipeline_output_folder)

    einverted_command = "einverted -sequence {} -gap 12 -threshold 50 -match 3 -mismatch -4 -outseq {} -outfile {}".format(assembly_fasta_path, einverted_output_fasta_path, einverted_alignments_path)
    gpf.run_system_command(einverted_command)
    gpf.run_system_command("convert_einverted_output_to_gff.py {} {} > {}".format(einverted_alignments_path, assembly_fasta_path, einverted_gff_path))

    gpf.run_system_command("gff_features_to_bedgraph.py {} inverted_repeat einverted_inverted_repeat --chunk_size {} > {}".format(einverted_gff_path, str(chunk_size), einverted_bedgraph_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder path for output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder path for main output files of the pipeline (where bedgraph files will be saved)")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.assembly_fasta_path, args.out_folder, args.pipeline_output_folder, args.chunk_size)

