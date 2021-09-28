#!/usr/bin/env python3
"""
Script for running OrthoMCL (including Diamond blastp for OrthoMCL)
"""

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse


def get_orthomcl_output_file_path(orthomcl_error_stream_path):
    """
    Function for getting OrthoMCL output file path from OrthoMCL's error stream file
    """
    orthomcl_error_stream_data = gpf.l(orthomcl_error_stream_path)
    orthomcl_file_path = None
    for line in orthomcl_error_stream_data:
        if "Final ORTHOMCL Result: " in line:
            orthomcl_file_path = gpf.spl(line, "Result: ", "generated")
            print(orthomcl_file_path)
            break
    if orthomcl_file_path is None:
        sys.stderr.write("Failed to find OrthoMCL output file path in the OrthoMCL error stream file (" + orthomcl_error_stream_path + ")\n")
        sys.exit(1)
    return orthomcl_file_path


def fix_characters_in_concat_fasta(concat_proteomes_temp_fasta_path, concat_proteomes_fasta_path):
    """
    Takes an input proteome file and outputs a proteome file where characters that break Diamond (#+.) have been replaced with an asterisk
    """
    with open(concat_proteomes_fasta_path, "w") as f:
        proteins_temp_data = gpf.ll(concat_proteomes_temp_fasta_path)
        for line in proteins_temp_data:
            if line.startswith(">") == False:
                line = line.replace("#", "*")
                line = line.replace("+", "*")
                line = line.replace(".", "*")
            f.write(line + "\n")
    os.remove(concat_proteomes_temp_fasta_path)


def main(proteomes_folder, csv_for_gg_file, orthomcl_folder, orthomcl_error_stream_path, threads, memory_limit, diamond_nonsensitive):
    gpf.run_system_command("mkdir -p " + orthomcl_folder)
    concat_proteomes_temp_fasta_path = orthomcl_folder + "/concat_proteomes_temp.faa"
    concat_proteomes_fasta_path = orthomcl_folder + "/concat_proteomes.faa"
    gpf.run_system_command("cat {}/*.f* > {}".format(proteomes_folder, concat_proteomes_temp_fasta_path))

    fix_characters_in_concat_fasta(concat_proteomes_temp_fasta_path, concat_proteomes_fasta_path)

    diamond_makedb_command = "diamond makedb --in {} -d {}".format(concat_proteomes_fasta_path, concat_proteomes_fasta_path)
    gpf.run_system_command(diamond_makedb_command)

    diamond_output_file = orthomcl_folder + "/diamond_blastp_output.m8" 

    diamond_sensitive_flag = " --sensitive "
    if diamond_nonsensitive == True:
        diamond_sensitive_flag = ""

    diamond_command = "diamond blastp -q {} -d {} -o {}{}--memory-limit {} --evalue 1e-5 --outfmt 6 --threads {}".format(concat_proteomes_fasta_path, concat_proteomes_fasta_path, diamond_output_file, diamond_sensitive_flag, str(memory_limit), str(threads))
    gpf.run_system_command(diamond_command)

    gg_file_path = orthomcl_folder + "/gg_file.txt"
    gg_file_command = "generate_orthomcl_gg_file_from_fasta.py {} {} > {}".format(csv_for_gg_file, proteomes_folder, gg_file_path)
    gpf.run_system_command(gg_file_command)

    orthomcl_command = "orthomcl.pl --mode 3 --blast_file {} --gg_file {} 2> {}".format(diamond_output_file, gg_file_path, orthomcl_error_stream_path)
    gpf.run_system_command(orthomcl_command)

    orthomcl_file_path = get_orthomcl_output_file_path(orthomcl_error_stream_path)

    move_orthomcl_file_command = "mv {} {}".format(orthomcl_file_path, orthomcl_folder)
    gpf.run_system_command(move_orthomcl_file_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proteomes_folder", type=str, help="Path to the protein FASTA files for OrthoMCL")
    parser.add_argument("csv_for_gg_file", type=str, help="Path to a CSV file that contains the information needed to generate a gg_file. First column: species identifiers (short names). Second column: names of the corresponding protein FASTA files for each species (without folder path)")
    parser.add_argument("orthomcl_folder", type=str, help="Path for the working folder for this script")
    parser.add_argument("orthomcl_error_stream_path", type=str, help="Path to the error stream file of OrthoMCL (this is needed to find the folder where OrthoMCL writes its files)")
    parser.add_argument("--threads", type=int, help="Number of threads for Diamond (default: 12)", default=12)
    parser.add_argument("--memory_limit", type=int, help="Memory limit for Diamond (in gigabytes, default: 5)", default=5)
    parser.add_argument("--diamond_nonsensitive", dest="diamond_nonsensitive", action="store_true", help="Optional: run Diamond blastp without the --sensitive flag to increase its running speed")
    args = parser.parse_args()
    main(args.proteomes_folder, args.csv_for_gg_file, args.orthomcl_folder, args.orthomcl_error_stream_path, args.threads, args.memory_limit, args.diamond_nonsensitive)