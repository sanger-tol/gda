#!/usr/bin/env python3
"""
Script for running tRNAscan-SE for detecting tRNAs
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
import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf


def main(pipeline_run_folder, assembly_fasta_path, threads):
    gene_annotation_folder = pipeline_run_folder + "/gene_annotation"
    gpf.run_system_command("mkdir -p " + gene_annotation_folder)
    trnascan_folder = gene_annotation_folder + "/trnascan"
    gpf.run_system_command("mkdir -p " + trnascan_folder)
    fasta_filename = assembly_fasta_path.split("/")[-1]
    fasta_basename = fasta_filename.split(".")[0]
    trnascan_bed_path = trnascan_folder + "/" + fasta_basename + "_tRNAscan_tRNAs.bed"
    trnascan_gff_path = trnascan_folder + "/" + fasta_basename + "_tRNAscan_tRNAs.gff3"
    trnascan_txt_path = trnascan_folder + "/" + fasta_basename + "_tRNAscan_tRNAs.txt"
    trnascan_stdout_path = trnascan_folder + "/" + fasta_basename + "_tRNAscan_stdout.txt"

    trnascan_command = "tRNAscan-SE {} --thread {} -o {} -b {} > {}".format(assembly_fasta_path, threads, trnascan_txt_path, trnascan_bed_path, trnascan_stdout_path)
    gpf.run_system_command(trnascan_command)

    trnascan_bed_to_gff_command = "convert_trnascan_bed_file_to_gff.py {} {} > {}".format(assembly_fasta_path, trnascan_bed_path, trnascan_gff_path)
    gpf.run_system_command(trnascan_bed_to_gff_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pipeline_run_folder", type=str, help="Folder for running the pipeline")
    parser.add_argument("assembly_fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    args = parser.parse_args()
    main(args.pipeline_run_folder, args.assembly_fasta_path, args.threads)
