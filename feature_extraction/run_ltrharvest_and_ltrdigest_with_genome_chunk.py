#!/usr/bin/env python3
"""
Script for running LTRharvest and LTRdigest
"""

import general_purpose_functions as gpf
import argparse
import os
import sys


def main(fasta_path, ltrharvest_folder):
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
    
    ltrdigest_temp_folder_path = ltrharvest_folder + "/temp_files"
    gpf.run_system_command("mkdir -p " + ltrdigest_temp_folder_path)
    os.chdir(ltrdigest_temp_folder_path)
    os.environ["TMPDIR"] = "."

    gpf.run_system_command("mkdir -p " + ltrharvest_folder)
    os.chdir(ltrharvest_folder)
    gpf.run_system_command("cp {} {}".format(fasta_path, ltrharvest_folder_fasta_path))

    gpf.run_system_command("gt suffixerator -db {} -indexname {} -tis -suf -lcp -des -ssp -sds -dna > {}".format(ltrharvest_folder_fasta_path, ltrharvest_folder_fasta_path, gt_suffixerator_stdout_path))
    gpf.run_system_command("gt ltrharvest -seqids -index {} -gff3 {} > {}".format(ltrharvest_folder_fasta_path, ltrharvest_gff_path, ltrharvest_stdout_path))
    if os.path.getsize(ltrharvest_gff_path) > 0:
        gpf.run_system_command("gt gff3 -sort {} > {}".format(ltrharvest_gff_path, ltrharvest_sorted_gff_path))
        gpf.run_system_command("gt ltrdigest -hmms {}/*.hmm -aaout -outfileprefix {} -seqfile {} -matchdescstart < {} > {}".format(hmms_folder, fasta_basename, ltrharvest_folder_fasta_path, ltrharvest_sorted_gff_path, ltrdigest_gff_path))
        gpf.run_system_command("gt select -rule_files {} -- < {} > {}".format(lua_file_path, ltrdigest_gff_path, filtered_ltrdigest_gff_path))
    else:
        sys.stderr.write("LTRharvest did not find any putative LTRs in input file {}\n".format(ltrharvest_folder_fasta_path))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("ltrharvest_folder", type=str, help="Path to folder for LTRharvest and LTRdigest output files")
    args = parser.parse_args()
    main(args.fasta_path, args.ltrharvest_folder)