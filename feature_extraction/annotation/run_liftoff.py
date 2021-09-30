#!/usr/bin/env python3
"""
Script for running Liftoff to transfer gene annotations
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


def gff_to_augustus_hints_file(gff_path, hints_file_path):
    """
    Filters a GFF to keep only exon and CDS features and reformats it
        so that it can be used as a hints file for Augustus
    """
    gff_data = gpf.ll(gff_path)
    with open(hints_file_path, "w") as hints_file:
        for line in gff_data:
            if line.startswith("#") == False:
                split_line = line.split()
                if len(split_line) > 2:
                    feature_type = split_line[2]
                    if feature_type == "exon" or feature_type == "CDS":
                        line += ";source=M"
                        hints_file.write(line + "\n")


def main(target_assembly_path, reference_assembly_path, reference_gff_path, output_gff_path, output_augustus_hints_file_path, threads):
    liftoff_command = "liftoff {} {} -g {} -p {} -o {}".format(target_assembly_path, reference_assembly_path, reference_gff_path, threads, output_gff_path)
    gpf.run_system_command(liftoff_command)
    gff_to_augustus_hints_file(output_gff_path, output_augustus_hints_file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_assembly_path", type=str, help="Path to the assembly that will be annotated")
    parser.add_argument("reference_assembly_path", type=str, help="Path to the reference assembly")
    parser.add_argument("reference_gff_path", type=str, help="Path to the reference GFF file")
    parser.add_argument("output_gff_path", type=str, help="Path for output GFF file")
    parser.add_argument("output_augustus_hints_file_path", type=str, help="Path for output Augustus hints file")
    parser.add_argument("--threads", type=int, help="Number of CPU threads (default: 1)", default=1)
    args = parser.parse_args()
    main(args.target_assembly_path, args.reference_assembly_path, args.reference_gff_path, args.output_gff_path, args.output_augustus_hints_file_path, args.threads)
