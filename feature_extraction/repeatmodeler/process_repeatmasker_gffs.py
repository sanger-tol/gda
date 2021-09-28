#!/usr/bin/env python3
"""
Script for running scripts that process RepeatMasker gff files
"""

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse

def main(assembly_title, fasta_path, repeatmasker_gff_path, pipeline_output_folder, chunk_size):
    repeatmasker_gff_basename = repeatmasker_gff_path.split("/")[-1]
    repeatmasker_gff_basename = repeatmasker_gff_basename.split(".fa.out.gff")[0]
    repeatmasker_gff_basename = repeatmasker_gff_basename.split(".fasta.out.gff")[0]

    complex_repeats_gff_path = pipeline_output_folder + "/" + repeatmasker_gff_basename + "_complex_repeats.gff"
    simple_repeats_gff_path = pipeline_output_folder + "/" + repeatmasker_gff_basename + "_simple_repeats.gff"
    simple_repeats_condensed_gff_path = pipeline_output_folder + "/" + repeatmasker_gff_basename + "_simple_repeats_condensed.gff"
    simple_repeats_csv_path = pipeline_output_folder + "/" + repeatmasker_gff_basename + "_simple_repeats.csv"

    gpf.run_system_command("mkdir -p " + pipeline_output_folder)

    os_command = "reformat_repeatmasker_gff.py " + fasta_path + " " + repeatmasker_gff_path + " complex " + complex_repeats_gff_path
    gpf.run_system_command(os_command)

    os_command = "reformat_repeatmasker_gff.py " + fasta_path + " " + repeatmasker_gff_path + " simple " + simple_repeats_gff_path
    gpf.run_system_command(os_command)

    os_command = "condense_simple_repeat_sequences.py " + simple_repeats_gff_path + " > " + simple_repeats_condensed_gff_path
    gpf.run_system_command(os_command)

    os_command = "repeatmasker_simple_repeat_frequencies.py " + simple_repeats_condensed_gff_path + " " + simple_repeats_csv_path
    gpf.run_system_command(os_command)

    simple_repeats_bedgraph_folder = pipeline_output_folder + "/simple_repeats_bedgraph"
    complex_repeats_bedgraph_folder = pipeline_output_folder + "/complex_repeats_bedgraph"

    os_command = "repeatmasker_gff_to_bedgraph.py {} {} {} {} --chunk_size {} --assembly_fasta_path {}".format(simple_repeats_condensed_gff_path, simple_repeats_bedgraph_folder, assembly_title, "simple_repeats", chunk_size, fasta_path)
    gpf.run_system_command(os_command)

    os_command = "repeatmasker_gff_to_bedgraph.py {} {} {} {} --chunk_size {} --assembly_fasta_path {}".format(complex_repeats_gff_path, complex_repeats_bedgraph_folder, assembly_title, "complex_repeats", chunk_size, fasta_path)
    gpf.run_system_command(os_command)

    fasta_basename = fasta_path.split("/")[-1]
    fasta_basename = fasta_basename.split(".")[0]
    simple_repeats_sum_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_simple_repeats_sum.bedgraph"
    complex_repeats_sum_bedgraph_path = pipeline_output_folder + "/" + fasta_basename + "_complex_repeats_sum.bedgraph"
    os_command = "sum_simple_or_complex_repeat_tracks.py {} {} {} {} --chunk_size {} > {}".format(simple_repeats_bedgraph_folder, fasta_path, assembly_title, "simple_repeats", chunk_size, simple_repeats_sum_bedgraph_path)
    gpf.run_system_command(os_command)
    os_command = "sum_simple_or_complex_repeat_tracks.py {} {} {} {} --chunk_size {} > {}".format(complex_repeats_bedgraph_folder, fasta_path, assembly_title, "complex_repeats", chunk_size, complex_repeats_sum_bedgraph_path)
    gpf.run_system_command(os_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_title", type=str, help="Name of the assembly")
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("repeatmasker_gff_path", type=str, help="Path to the GFF output from RepeatMasker")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for output files")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.assembly_title, args.fasta_path, args.repeatmasker_gff_path, args.pipeline_output_folder, str(args.chunk_size))

