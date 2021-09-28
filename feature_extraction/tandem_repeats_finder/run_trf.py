#!/usr/bin/env python3
"""
Script for running Tandem Repeats Finder as a part of genome decomposition
"""

import sys
import os
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import argparse
import general_purpose_functions as gpf


def mask_Ns_in_fasta(in_path, out_path):
    """
    Input: 1) path to input FASTA file that may contain Ns, 2) path for output FASTA file, where Ns will be replaced with Xs
    Output: new FASTA file where Ns have been replaced with Xs
    """
    out_list = list()
    fasta_data = gpf.ll(in_path)
    for line in fasta_data:
        if line.startswith(">") == False:
            line = line.upper()
            line = line.replace("N", "X")
        out_list.append(line)
    gpf.export_list_as_line_break_separated_file(out_list, out_path)


def main(in_path, out_folder, pipeline_output_folder, chunk_size):
    fasta_filename = in_path.split("/")[-1]
    fasta_basename = fasta_filename.split(".")[0]

    fasta_with_masked_Ns_path = out_folder + "/" + fasta_basename + "_masked_N.fa"
    mask_Ns_in_fasta(in_path, fasta_with_masked_Ns_path)
    gpf.run_system_command("mkdir -p " + out_folder)
    gpf.run_system_command("mkdir -p " + pipeline_output_folder)
    os.chdir(out_folder)
    trf_command = "trf " + fasta_with_masked_Ns_path + " 2 1000 1000 80 10 25 1000 -m -h -ngs > /dev/null" 
    gpf.run_system_command(trf_command)
    trf_out_file_path = fasta_with_masked_Ns_path + ".2.1000.1000.80.10.25.1000.mask"

    repeat_density_out_path = fasta_basename + "_tandem_repeat_density.txt"
    repeat_density_command = "trf_repeat_density_sliding_window.py " + trf_out_file_path + " --chunk_size " + str(chunk_size) + " --zero_based > " + repeat_density_out_path
    gpf.run_system_command(repeat_density_command)

    repeat_density_gff_out_path = fasta_basename + "_tandem_repeat_density.gff"
    repeat_density_bed_out_path = fasta_basename + "_tandem_repeat_density.bed"
    bedgraph_out_path = pipeline_output_folder + "/" + fasta_basename + "_tandem_repeat_density.bedgraph"
    repeat_density_gff_command = "trf_repeat_density_to_gff.py " + repeat_density_out_path + " GFF --assembly_path " + in_path + " > " + repeat_density_gff_out_path
    repeat_density_bed_command = "trf_repeat_density_to_gff.py " + repeat_density_out_path + " BED > " + repeat_density_bed_out_path
    repeat_density_bedgraph_command = "trf_repeat_density_to_bedgraph.py " + repeat_density_out_path + " > " + bedgraph_out_path
    gpf.run_system_command(repeat_density_gff_command)
    gpf.run_system_command(repeat_density_bed_command)
    gpf.run_system_command(repeat_density_bedgraph_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to input FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder for output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph file")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk length (bp) for the estimation of repeat density using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.in_path, args.out_folder, args.pipeline_output_folder, args.chunk_size)
    
