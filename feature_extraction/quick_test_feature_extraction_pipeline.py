#!/usr/bin/env python3
"""
Script for running a quick test of the genomic feature extraction pipeline of GDA with P. falciparum chromosome 1 as the input.
This script runs all the mandatory parts of the pipeline but skips most of the optional parts
"""

import pandas as pd
import os
import argparse
import general_purpose_functions as gpf
import sys


def compare_test_df_with_ref_df(test_df_path, ref_df_path):
    """
    Function for comparing test run results with precomputed refrerence results. Returns true if the results match, false if they don't
    """
    test_df = pd.read_csv(test_df_path, sep="\t")
    ref_df = pd.read_csv(ref_df_path, sep="\t")
    test_pass_flag = True
    ref_df_colnames = list(ref_df.columns)
    for ref_df_colname in ref_df_colnames:
        ref_col = ref_df[ref_df_colname]
        test_col = test_df[ref_df_colname]
        if ref_col.equals(test_col) == False:
            test_pass_flag = False
            break
    return test_pass_flag


def main(test_run_folder, singularity_image_path):
    test_run_folder = os.path.abspath(test_run_folder)
    if singularity_image_path != "":
        singularity_image_path = os.path.abspath(singularity_image_path)

    if singularity_image_path != "" and os.path.isfile(singularity_image_path) == False:
        sys.stderr.write("Singularity image file was not found at {}\n".format(singularity_image_path))
        sys.exit(1) 

    gpf.run_system_command("mkdir -p {}".format(test_run_folder))
    if len(os.listdir(test_run_folder)) > 0:
        sys.stderr.write("Test run folder ({}) exists and is not empty\n".format(test_run_folder))
        sys.exit(1)

    current_script_folder = os.path.dirname(os.path.realpath(__file__))
    split_current_script_folder = current_script_folder.split("/")
    test_data_folder = "/".join(split_current_script_folder[0:len(split_current_script_folder) - 1]) + "/test_data"
    test_fasta_path = test_data_folder + "/Pf3D7_chr1.fa"
    ref_df_path = test_data_folder + "/Pf3D7_chr1_merged_bedgraph.tsv"
    print("test_fasta_path", test_fasta_path)

    if os.path.isfile(test_fasta_path) == False:
        sys.stderr.write("Test FASTA file was not found at {}\n".format(test_fasta_path))
        sys.exit(1)

    if os.path.isfile(ref_df_path) == False:
        sys.stderr.write("Precomputed reference TSV file was not found at {}\n".format(ref_df_path))
        sys.exit(1)

    test_pass_flag = True
    singularity_string = ""
    if singularity_image_path != "":
        singularity_string = "--singularity_image_path {} ".format(singularity_image_path)
    test_run_command = "gda extract_genomic_features {} --threads 1 {}--telomeric_seq_preset apicomplexan --pipeline_run_folder {} --run_gene_annotation_pipeline --augustus_species pfalciparum".format(test_fasta_path, singularity_string, test_run_folder)
    test_run_exit_code = os.system(test_run_command)
    if test_run_exit_code != 0:
        test_pass_flag = False
    else:
        test_df_path = test_run_folder + "/merged_bedgraph_table/Pf3D7_chr1_merged_bedgraph.tsv"
        if compare_test_df_with_ref_df(test_df_path, ref_df_path) == False:
            test_pass_flag = False
    if test_pass_flag == True:
        print("Successfully completed the test run")
    else:
        print("The test run failed")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_run_folder", type=str, help="Folder for executing the test run")
    parser.add_argument("--singularity_image_path", type=str, help="Path to GDA Singularity image (needed if GDA dependencies are not already in path)", default="")
    args = parser.parse_args()
    main(args.test_run_folder, args.singularity_image_path)