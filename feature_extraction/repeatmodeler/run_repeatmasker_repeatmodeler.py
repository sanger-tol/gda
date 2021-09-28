#!/usr/bin/env python3
"""
Script for running RepeatMasker and RepeatModeler as a part of genome decomposition
"""

import os
import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse


def check_fasta_header_lengths(fasta_path, max_header_len=50):
    """"
    Function for checking if the FASTA headers are short enough for RepeatMasker, as RepeatMasker requires the headers to be at most 50 characters long
    """
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    for fasta_tuple in fasta_data:
        header = str(fasta_tuple[0])
        if len(header) > max_header_len:
            sys.stderr.write("Error: FASTA header is too long (the maximum length is " + str(max_header_len) + " characters)\n")
            sys.stderr.write("FASTA header: " + header + "\n")
            sys.exit(1)


def check_if_repeat_families_detected(repeatmodeler_stdout_path):
    """
    Checks the STDOUT of RepeatModeler for the error message that appears when no repeat families are detected
    """
    repeat_families_detected_flag = True
    repeatmodeler_stdout_data = gpf.l(repeatmodeler_stdout_path)
    for line in repeatmodeler_stdout_data:
        if line == "No families identified.  Perhaps the database is too small":
            repeat_families_detected_flag = False
            break
    return repeat_families_detected_flag


def main(assembly_title, fasta_path, repeatmodeler_folder, pipeline_output_folder, repeatmodeler_recovery_folder, threads, chunk_size):
    check_fasta_header_lengths(fasta_path)

    repeatmodeler_db_folder = repeatmodeler_folder + "/repeatmodeler_db"
    fasta_filename = fasta_path.split("/")[-1]
    repeatmodeler_db_path = repeatmodeler_db_folder + "/" + fasta_filename
    repeatmasker_folder = repeatmodeler_folder + "/repeatmasker"
    repeatmasker_fasta_path = repeatmasker_folder + "/" + fasta_filename
    repeatmodeler_repeatmasker_stdout_path = repeatmasker_folder + "/repeatmodeler_repeatmasker_stdout.txt"

    gpf.run_system_command("mkdir -p " + repeatmodeler_folder)
    gpf.run_system_command("mkdir -p " + repeatmodeler_db_folder)
    gpf.run_system_command("mkdir -p " + repeatmasker_folder)
    os.chdir(repeatmodeler_db_folder)

    if repeatmodeler_recovery_folder == "NA":
        builddb_command = "BuildDatabase -name {} -engine ncbi {} > {}".format(fasta_filename, fasta_path, repeatmodeler_repeatmasker_stdout_path)
        gpf.run_system_command(builddb_command)
    os.chdir(repeatmodeler_folder)

    repeatmodeler_command = "RepeatModeler -engine ncbi -pa 3 -srand 1 -database " + repeatmodeler_db_path
    if repeatmodeler_recovery_folder != "NA":
        repeatmodeler_command += " -recoverDir " + repeatmodeler_recovery_folder
    repeatmodeler_command += " >> {}".format(repeatmodeler_repeatmasker_stdout_path)
    gpf.run_system_command(repeatmodeler_command)

    repeat_families_detected_flag = check_if_repeat_families_detected(repeatmodeler_repeatmasker_stdout_path)
    if repeat_families_detected_flag == True:

        repeatmodeler_subfolders = [f.name for f in os.scandir(repeatmodeler_folder) if f.is_dir()]
        repeatmodeler_output_folders = [n for n in repeatmodeler_subfolders if n.startswith("RM_")]
        repeatmodeler_output_folder = None
        repeatmodeler_output_folders_len = len(repeatmodeler_output_folders)
        if repeatmodeler_output_folders_len == 1:
            repeatmodeler_output_folder = repeatmodeler_output_folders[0]
        elif repeatmodeler_output_folders_len == 0:
            sys.stderr.write("Error: could not find RepeatModeler output folder\n")
            sys.exit(1)
        else:
            sys.stderr.write("Error: found more than 1 folder when looking for RepeatModeler output folder\n")
            sys.exit(1)

        repeatmodeler_consensi_path = repeatmodeler_folder + "/" + repeatmodeler_output_folder + "/consensi.fa"
        if os.path.isfile(repeatmodeler_consensi_path) == False:
            sys.stderr.write("Error: could not find consensi.fa output file of RepeatModeler\n")
            sys.exit(1)

        os.chdir(repeatmasker_folder)

        gpf.run_system_command("cp " + fasta_path + " " + repeatmasker_fasta_path)

        repeatmasker_command = "RepeatMasker -pa {} -gff -e rmblast -lib {} {} >> {}".format(str(threads), repeatmodeler_consensi_path, repeatmasker_fasta_path, repeatmodeler_repeatmasker_stdout_path)
        gpf.run_system_command(repeatmasker_command)

        repeatmasker_gff_path = repeatmasker_folder + "/" + fasta_filename + ".out.gff"
        if os.path.isfile(repeatmasker_gff_path) == False:
            sys.stderr.write("Error: could not find RepeatMasker output GFF file\n")
            sys.exit(1)
        else:
            print("repeatmasker_gff_path:", repeatmasker_gff_path)

        gff_processing_command = "process_repeatmasker_gffs.py {} {} {} {} --chunk_size {} >> {}".format(assembly_title, fasta_path, repeatmasker_gff_path, pipeline_output_folder, chunk_size, repeatmodeler_repeatmasker_stdout_path)
        gpf.run_system_command(gff_processing_command)
    
    else:
        sys.stderr.write("No repeat families were detected by RepeatModeler. RepeatModeler exited with this message:\n")
        sys.stderr.write("'No families identified. Perhaps the database is too small or contains overly fragmented sequences.'\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_title", type=str, help="Name of the assembly")
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("repeatmodeler_folder", type=str, help="Folder for running RepeatModeler")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for output files")
    parser.add_argument("--repeatmodeler_recovery_folder", type=str, help="Optional: path to a previous unfinished run of RepeatModeler (for recovering a crashed run), default=NA", default="NA")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.assembly_title, args.fasta_path, args.repeatmodeler_folder, args.pipeline_output_folder, args.repeatmodeler_recovery_folder, args.threads, str(args.chunk_size))