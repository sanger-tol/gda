#!/usr/bin/env python3
"""
Script for splitting a genome assembly into chunks with fixed size and running LTRharvest + LTRdigest with the chunks
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
import os
import argparse
import sys
import shutil


def split_assembly_to_chunks(run_folder, input_fasta_path, fasta_chunk_size):
    """
    Splits the input assembly FASTA file into smaller chunk files with fixed size
    """
    assembly_fasta_path = run_folder + "/assembly.fasta"
    gpf.run_system_command("shorten_fasta_headers.py {} > {}".format(input_fasta_path, assembly_fasta_path))
    gpf.run_system_command("samtools faidx {}".format(assembly_fasta_path))

    seq_len_file_path = run_folder + "/seq_len_file.txt"
    gpf.run_system_command("get_fasta_sequence_lengths.py --tab {} > {}".format(assembly_fasta_path, seq_len_file_path))
    windows_bed_file_path = run_folder + "/fasta_windows.bed"

    gpf.run_system_command("bedtools makewindows -g {} -w {} -i srcwinnum > {}".format(seq_len_file_path, fasta_chunk_size, windows_bed_file_path))
    fasta_windows_path = run_folder + "/fasta_windows.fa"

    gpf.run_system_command("bedtools getfasta -fi {} -fo {} -bed {}".format(assembly_fasta_path, fasta_windows_path, windows_bed_file_path))

    split_fasta_folder = run_folder + "/split_fasta"
    gpf.run_system_command("mkdir -p " + split_fasta_folder)
    gpf.run_system_command("split_multifasta_with_fixed_nr_of_sequences.py {} {} 1".format(fasta_windows_path, split_fasta_folder))



def run_ltrharvest_ltrdigest_with_chunks(run_folder, threads):
    """
    Runs LTRharvest and LTRdigest with genome assembly chunks
    """
    split_fasta_folder = run_folder + "/split_fasta"
    split_fasta_files = sorted(gpf.get_file_paths(split_fasta_folder, "fa"))

    ltrharvest_ltrdigest_folder = run_folder + "/ltrharvest_ltrdigest_chunks"
    gpf.run_system_command("mkdir -p " + ltrharvest_ltrdigest_folder)

    parallel_args = list()
    for split_fasta_file in split_fasta_files:
        split_fasta_basename = split_fasta_file.split("/")[-1]
        split_fasta_basename = split_fasta_basename[0:len(split_fasta_basename) -3]
        split_fasta_folder = ltrharvest_ltrdigest_folder + "/" + split_fasta_basename
        gpf.run_system_command("mkdir -p " + split_fasta_folder)
        parallel_args.append(split_fasta_file)
        parallel_args.append(split_fasta_folder)
    parallel_args_file_path = ltrharvest_ltrdigest_folder + "/parallel_args.txt"
    gpf.export_list_as_line_break_separated_file(parallel_args, parallel_args_file_path)

    parallel_command = "parallel -j {} --max-args 2 -a {} run_ltrharvest_and_ltrdigest_with_genome_chunk.py".format(threads, parallel_args_file_path)
    gpf.run_system_command(parallel_command)


def get_ltrdigest_gff_paths(ltrdigest_chunks_folder):
    """
    Finds the file paths of GFF files that contain filtered LTRdigest matches for genome chunks
    """
    ltrdigest_subfolders = [f.path for f in os.scandir(ltrdigest_chunks_folder) if f.is_dir()]
    ltrdigest_filtered_gffs = list()
    for ltrdigest_subfolder in ltrdigest_subfolders:
        gff_files = gpf.get_file_paths(ltrdigest_subfolder, "gff")
        gff_files = [n for n in gff_files if n.endswith("_ltrdigest_filtered.gff")]
        assert len(gff_files) < 2
        if len(gff_files) == 1:
            ltrdigest_filtered_gff_path = gff_files[0]
            ltrdigest_filtered_gffs.append(ltrdigest_filtered_gff_path)
    return ltrdigest_filtered_gffs


def merge_ltrdigest_gffs(ltrdigest_filtered_gffs, output_gff_path):
    """
    Reads the LTRdigest output GFF files of genome chunks and merges them into one GFF file for the full genome
    """
    with open(output_gff_path, "w") as f:
        f.write("##gff-version 3\n")

        for ltrdigest_filtered_gff_path in ltrdigest_filtered_gffs:
            if os.path.getsize(ltrdigest_filtered_gff_path) > 0:
                gff_data = gpf.l(ltrdigest_filtered_gff_path)
                assert gff_data[0] == "##gff-version 3"
                seq_region_line = gff_data[1]
                split_seq_region_line = seq_region_line.split()
                assert len(split_seq_region_line) == 4
                region_coords = split_seq_region_line[1].split(":")
                assert len(region_coords) == 2
                scaff_name = region_coords[0]
                region_offset = int(region_coords[1].split("-")[0])

                for line in gff_data[2:len(gff_data)]:
                    if line.startswith("#") == False:
                        split_line = line.split()
                        assert len(split_line) == 9
                        split_line[0] = scaff_name
                        start_coord = int(split_line[3]) + region_offset
                        end_coord = int(split_line[4]) + region_offset
                        split_line[3] = str(start_coord)
                        split_line[4] = str(end_coord)
                        edited_line = "\t".join(split_line)
                        f.write(edited_line + "\n")
            else:
                sys.stderr.write("Filtered LTRdigest GFF file ({}) is empty\n".format(ltrdigest_filtered_gff_path))


def main(input_fasta_path, run_folder, pipeline_output_folder, bedgraph_chunk_size, fasta_chunk_size, threads):
    input_fasta_path = os.path.abspath(input_fasta_path)
    run_folder = os.path.abspath(run_folder)
    pipeline_output_folder = os.path.abspath(pipeline_output_folder)

    gpf.run_system_command("mkdir -p " + run_folder)
    if len(os.listdir(run_folder)) > 0:
        sys.stderr.write("Warning: the directory specified as run_folder ({}) already exists and is not empty\n".format(run_folder))
        #sys.exit(1)

    gpf.run_system_command("mkdir -p " + pipeline_output_folder)
    os.chdir(run_folder)

    split_assembly_to_chunks(run_folder, input_fasta_path, fasta_chunk_size)
    run_ltrharvest_ltrdigest_with_chunks(run_folder, threads)
    ltrdigest_chunks_folder = run_folder + "/ltrharvest_ltrdigest_chunks"
    #gpf.run_system_command("mkdir -p " + ltrdigest_chunks_folder)

    fasta_basename_with_extension = input_fasta_path.split("/")[-1]
    fasta_basename = fasta_basename_with_extension.split(".")[0]
    merged_gff_path = run_folder + "/ltrdigest_chunks_merged.gff"

    ltrdigest_filtered_gffs = get_ltrdigest_gff_paths(ltrdigest_chunks_folder)
    merge_ltrdigest_gffs(ltrdigest_filtered_gffs, merged_gff_path)
    if os.path.exists(merged_gff_path) == False:
        sys.stderr.write("Failed to generate a merged GFF file from LTRdigest output files\n")
        sys.exit(1)

    ltr_retrotransposon_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_ltrdigest_retrotransposons.bedgraph"
    ltr_protein_match_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_ltrdigest_protein_matches.bedgraph"

    gpf.run_system_command("gff_features_to_bedgraph.py {} LTR_retrotransposon LTRdigest_LTR_retrotransposon --chunk_size {} --assembly_fasta_path {} > {}".format(merged_gff_path, bedgraph_chunk_size, input_fasta_path, ltr_retrotransposon_bedgraph_path))
    gpf.run_system_command("gff_features_to_bedgraph.py {} protein_match LTRdigest_protein_match --chunk_size {} --assembly_fasta_path {} > {}".format(merged_gff_path, bedgraph_chunk_size, input_fasta_path, ltr_protein_match_bedgraph_path))

    split_fasta_folder = run_folder + "/split_fasta"
    shutil.rmtree(split_fasta_folder)
    shutil.rmtree(ltrdigest_chunks_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("run_folder", type=str, help="Path to folder for intermediate files for LTRharvest and LTRdigest on genome chunks")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph file")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--fasta_chunk_size", type=int, help="Length of sequence chunks (bp) into which the input FASTA file will be broken before running LTRharvest (default: 500000)", default=500000)
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    args = parser.parse_args()
    main(args.input_fasta_path, args.run_folder, args.pipeline_output_folder, args.chunk_size, args.fasta_chunk_size, args.threads)
