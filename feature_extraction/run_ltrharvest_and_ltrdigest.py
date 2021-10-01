#!/usr/bin/env python3
"""
Script for running LTRharvest and LTRdigest
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


def main(fasta_path, ltrharvest_folder, pipeline_output_folder, chunk_size):
    script_folder = os.path.dirname(os.path.abspath(__file__))
    lua_file_path = script_folder + "/third_party_files/filter_protein_match.lua"
    hmms_folder = script_folder + "/third_party_files/ltr_hmm"

    fasta_basename_with_extension = fasta_path.split("/")[-1]
    ltrharvest_folder_fasta_path = ltrharvest_folder + "/" + fasta_basename_with_extension
    fasta_basename = fasta_basename_with_extension.split(".")[0]

    ltrharvest_gff_path = ltrharvest_folder + "/" + fasta_basename + "_ltrharvest.gff"
    ltrharvest_sorted_gff_path = ltrharvest_folder + "/" + fasta_basename + "_ltrharvest_sorted.gff"
    ltrharvest_stdout_path = ltrharvest_folder + "/" + fasta_basename + "_ltrharvest_stdout.txt"
    gt_suffixerator_stdout_path = ltrharvest_folder + "/" + fasta_basename + "_gt_suffixerator_stdout.txt"
    ltrdigest_gff_path = ltrharvest_folder + "/" + fasta_basename + "_ltrdigest.gff"
    filtered_ltrdigest_gff_path = ltrharvest_folder + "/" + fasta_basename + "_ltrdigest_filtered.gff"
    ltr_retrotransposon_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_ltrdigest_retrotransposons.bedgraph"
    ltr_protein_match_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_ltrdigest_protein_matches.bedgraph"

    ltrdigest_temp_folder_path = ltrharvest_folder + "/temp_files"
    gpf.run_system_command("mkdir -p " + ltrdigest_temp_folder_path)
    os.chdir(ltrdigest_temp_folder_path)
    os.environ["TMPDIR"] = "."

    gpf.run_system_command("mkdir -p " + ltrharvest_folder)
    #os.chdir(ltrharvest_folder)
    gpf.run_system_command("cp {} {}".format(fasta_path, ltrharvest_folder_fasta_path))

    gpf.run_system_command("gt suffixerator -db {} -indexname {} -tis -suf -lcp -des -ssp -sds -dna > {}".format(ltrharvest_folder_fasta_path, ltrharvest_folder_fasta_path, gt_suffixerator_stdout_path))
    gpf.run_system_command("gt ltrharvest -seqids -index {} -gff3 {} > {}".format(ltrharvest_folder_fasta_path, ltrharvest_gff_path, ltrharvest_stdout_path))
    gpf.run_system_command("gt gff3 -sort {} > {}".format(ltrharvest_gff_path, ltrharvest_sorted_gff_path))
    #gpf.run_system_command("export TMPDIR=\"{}\"; gt ltrdigest -hmms {}/*.hmm -aaout -outfileprefix {} -seqfile {} -matchdescstart < {} > {}".format(ltrdigest_temp_folder_path, hmms_folder, fasta_basename, ltrharvest_folder_fasta_path, ltrharvest_sorted_gff_path, ltrdigest_gff_path))

    gpf.run_system_command("gt ltrdigest -hmms {}/*.hmm -aaout -outfileprefix {} -seqfile {} -matchdescstart < {} > {}".format(hmms_folder, fasta_basename, ltrharvest_folder_fasta_path, ltrharvest_sorted_gff_path, ltrdigest_gff_path))

    gpf.run_system_command("gt select -rule_files {} -- < {} > {}".format(lua_file_path, ltrdigest_gff_path, filtered_ltrdigest_gff_path))
    gpf.run_system_command("gff_features_to_bedgraph.py {} LTR_retrotransposon LTRdigest_LTR_retrotransposon --chunk_size {} > {}".format(filtered_ltrdigest_gff_path, str(chunk_size), ltr_retrotransposon_bedgraph_path))
    gpf.run_system_command("gff_features_to_bedgraph.py {} protein_match LTRdigest_protein_match --chunk_size {} > {}".format(filtered_ltrdigest_gff_path, str(chunk_size), ltr_protein_match_bedgraph_path))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("ltrharvest_folder", type=str, help="Path to folder for LTRharvest and LTRdigest output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph file")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.fasta_path, args.ltrharvest_folder, args.pipeline_output_folder, args.chunk_size)
