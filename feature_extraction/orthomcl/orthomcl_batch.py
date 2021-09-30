#!/usr/bin/env python3
"""
Script for running OrthoMCL as batch
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

import os
import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse
from glob import glob

def main(target_assembly_path, target_gff_path, orthomcl_folder, pipeline_output_folder, target_species_id, orthomcl_references_folder, threads, memory_limit, chunk_size, diamond_nonsensitive):
    fasta_basename_with_extension = target_assembly_path.split("/")[-1]
    fasta_basename = fasta_basename_with_extension.split(".")[0]

    gpf.run_system_command("mkdir -p " + orthomcl_folder)

    run_dir_references_folder = orthomcl_folder + "/references"
    copy_references_folder_command = "cp -ar {} {}".format(orthomcl_references_folder, run_dir_references_folder)
    gpf.run_system_command(copy_references_folder_command)
    target_proteome_path = orthomcl_folder + "/" + fasta_basename + "_proteome.faa"
    filtered_target_gff_path = orthomcl_folder + "/" + fasta_basename + "_filtered.gff3"
    gff_filtering_command = "remove_non_mrna_cds_features.py {} > {}".format(target_gff_path, filtered_target_gff_path)
    gpf.run_system_command(gff_filtering_command)
    proteome_extraction_command = "gff_to_transcripts_and_proteins.py {} {} {}".format(target_assembly_path, filtered_target_gff_path, orthomcl_folder)
    gpf.run_system_command(proteome_extraction_command)

    references_subfolders = glob("{}/*/".format(run_dir_references_folder))
    for references_subfolder in references_subfolders:
        references_subfolder = references_subfolder[0:len(references_subfolder) - 1]
        os.chdir(references_subfolder)
        link_reference_proteome_command = "ln -s {}".format(target_proteome_path)
        gpf.run_system_command(link_reference_proteome_command)
        table_for_gg_file_path = references_subfolder + "/table_for_gg_file.csv"
        target_proteome_filename = target_proteome_path.split("/")[-1]
        reference_set_id = references_subfolder.split("/")[-1]
        with open(table_for_gg_file_path, "a") as f:
            f.write("{},{}\n".format(target_species_id, target_proteome_filename))
        orthomcl_run_folder = references_subfolder + "/orthomcl_run"
        orthomcl_error_stream_path = orthomcl_run_folder + "/orthomcl_error_stream.txt"

        diamond_nonsensitive_flag = ""
        if diamond_nonsensitive == True:
            diamond_nonsensitive_flag = " --diamond_nonsensitive "

        orthomcl_run_command = "run_orthomcl.py {} {} {} {} --threads {} --memory_limit {}{}".format(references_subfolder, table_for_gg_file_path, orthomcl_run_folder, orthomcl_error_stream_path, threads, memory_limit, diamond_nonsensitive_flag)
        gpf.run_system_command(orthomcl_run_command)
        orthomcl_output_file_path = orthomcl_run_folder + "/all_orthomcl.out"
        orthomcl_conservation_file_path = orthomcl_run_folder + "/{}_orthomcl_conservation.csv".format(reference_set_id)
        orthomcl_conservation_command = "orthomcl_conservation.py {} {} {} {} --query_gff_feature {}".format(target_species_id, target_gff_path, orthomcl_output_file_path, orthomcl_conservation_file_path, "mRNA")
        gpf.run_system_command(orthomcl_conservation_command)
        orthomcl_to_bedgraph_command = "stats_per_gene_to_bedgraph.py {} {} OrthoMCL {} {} --reference_dataset_title {} --chunk_size {}".format(target_assembly_path, orthomcl_conservation_file_path, pipeline_output_folder, target_species_id, reference_set_id, chunk_size)
        gpf.run_system_command(orthomcl_to_bedgraph_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_assembly_path", type=str, help="Path to the target assembly")
    parser.add_argument("target_gff_path", type=str, help="Path to the GFF3 annotations file for the target assembly")
    parser.add_argument("orthomcl_folder", type=str, help="Path for the working folder for this script")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder path for main output files of the pipeline (where bedgraph files will be saved)")
    parser.add_argument("target_species_id", type=str, help="Species ID of the target species")
    parser.add_argument("orthomcl_references_folder", type=str, help="Path to the folder with the reference files for this script")
    parser.add_argument("--threads", type=int, help="Number of threads for Diamond (default: 1)", default=1)
    parser.add_argument("--memory_limit", type=int, help="Memory limit for Diamond (in gigabytes, default: 5)", default=5)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--diamond_nonsensitive", dest="diamond_nonsensitive", action="store_true", help="Optional: run Diamond blastp without the --sensitive flag to increase its running speed")
    args = parser.parse_args()
    main(args.target_assembly_path, args.target_gff_path, args.orthomcl_folder, args.pipeline_output_folder, args.target_species_id, args.orthomcl_references_folder, args.threads, args.memory_limit, args.chunk_size, args.diamond_nonsensitive)
