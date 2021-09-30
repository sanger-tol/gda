#!/usr/bin/env python3
"""
Script for extracting transcripts and protein sequences from a GFF3 and a genome assembly FASTA file
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
import os
import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf

def main(assembly_fasta_path, gff_path, out_folder):

    assembly_fasta_filename = assembly_fasta_path.split("/")[-1]
    assembly_fasta_basename = assembly_fasta_filename.split(".")[0]

    transcripts_path = out_folder + "/" + assembly_fasta_basename + "_transcripts.fa"
    proteins_temp_path = out_folder + "/" + assembly_fasta_basename + "_proteome_temp.faa"
    proteins_path = out_folder + "/" + assembly_fasta_basename + "_proteome.faa"

    gffread_command = "gffread --no-pseudo -w {} -g {} {}".format(transcripts_path, assembly_fasta_path, gff_path)
    gpf.run_system_command(gffread_command)

    gffread_command2 = "gffread --no-pseudo -y {} -g {} {}".format(proteins_temp_path, assembly_fasta_path, gff_path)
    gpf.run_system_command(gffread_command2)

    with open(proteins_path, "w") as f:
        proteins_temp_data = gpf.ll(proteins_temp_path)
        for line in proteins_temp_data:
            if line.startswith(">"):
                header = line[1:len(line)]
                header = "".join([c if (c.isalnum() or c=="." or c=="-") else "_" for c in header]) # Filters FASTA header to replace characters that are neither alphanumeric, dots or hyphens with underscores
                line = ">" + header
            else:
                if "." in line:
                    line = line.replace(".", "*")
            f.write(line + "\n")

    os.remove(proteins_temp_path)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_fasta_path", type=str, help="Path to input assembly FASTA file")
    parser.add_argument("gff_path", type=str, help="Path to gene annotations GFF file")
    parser.add_argument("out_folder", type=str, help="Path to folder for output files")
    args = parser.parse_args()
    main(args.assembly_fasta_path, args.gff_path, args.out_folder)


