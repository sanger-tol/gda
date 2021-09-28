#!/usr/bin/env python3
"""
Script for automatically running the gda_parameters.py and gda_clustering.py scripts after running the genomic feature extraction
"""

import argparse
import general_purpose_functions as gpf
import sys
import os.path


def get_selected_tsv_file_name(merged_tsv_folder):
    """
    Looks for TSV files in the specified folder and selects the one to be used for clustering, based on the file name
    """
    tsv_files = gpf.get_file_paths(merged_tsv_folder, "tsv")
    selected_tsv_file = None
    tsv_files_count = len(tsv_files)
    if tsv_files_count == 0:
        sys.stderr.write("No TSV files were found in the directory {}\n".format(merged_tsv_folder))
        sys.exit(1)
    elif tsv_files_count == 1:
        selected_tsv_file = tsv_files[0]
    elif tsv_files_count == 2:
        for tsv_file in tsv_files:
            if "_downsampled_" in tsv_file and tsv_file.endswith("_x.tsv"):
                selected_tsv_file = tsv_file
    else:
        sys.stderr.write("1 or 2 TSV files expected in the directory {} but {} was found\n".format(merged_tsv_folder, tsv_files_count))
        sys.exit(1)
    return selected_tsv_file


def main(merged_tsv_folder, out_folder):
    selected_tsv_file = get_selected_tsv_file_name(merged_tsv_folder)
    gda_parameters_command = "gda_parameters.py -d {} {}".format(out_folder, selected_tsv_file)
    gpf.run_system_command(gda_parameters_command)
    clustering_metrics_file_path = out_folder + "/parameter_selection/clustering_metrics.csv"
    if os.path.isfile(clustering_metrics_file_path):
        clustering_command = "gda_clustering.py -d {} -m {} {}".format(out_folder, clustering_metrics_file_path, selected_tsv_file)
        gpf.run_system_command(clustering_command)

    else:
        sys.stderr.write("Clustering metrics file ({}) was not found\n".format(clustering_metrics_file_path))
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("merged_tsv_folder", type=str, help="Path to the directory with the TSV file(s) produced from bedgraph files")
    parser.add_argument("out_folder", type=str, help="Path to the directory with clustering output files")
    args = parser.parse_args()
    main(args.merged_tsv_folder, args.out_folder)