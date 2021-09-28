#!/usr/bin/env python3
"""
Wrapper script for starting GDA Shiny app
"""

import argparse
import os
import sys


def main(in_folder, selected_species):
    in_folder = os.path.abspath(in_folder)
    if os.path.exists(in_folder) == False:
        sys.stderr.write("Folder {} was not found\n".format(in_folder))
        sys.exit(1)
    if os.path.isdir(in_folder) == False:
        sys.stderr.write("{} does not appear to be a directory\n".format(in_folder))
        sys.exit(1)
    in_folder_files = os.listdir(in_folder)
    expected_files = ("cluster_heatmap.csv", "feathist.json", "genomeprops.csv", "umap_clustering.csv", "clusters.gff", "feattable.csv")
    for expected_file in expected_files:
        if expected_file not in in_folder_files:
            sys.stderr.write("File {} was not found. The user-provided input directory ({}) does not appear to be a gda_out directory\n".format(expected_file, in_folder))
            sys.exit(1)
    current_script_dir = os.path.dirname(os.path.realpath(__file__))
    subfolders = [f.path for f in os.scandir(in_folder) if f.is_dir()]
    subfolders = [os.path.basename(os.path.normpath(n)) for n in subfolders]
    subfolders = [n for n in subfolders if n != "parameter_selection"]
    if selected_species == "":
        if len(subfolders) == 1:
            selected_species = subfolders[0]
        elif len(subfolders) > 1:
            sys.stderr.write("The gda_out folder appears to contain the clustering results of more than one species ({}) but no species was selected by the user with the --selected_species flag\n".format(", ".join(subfolders)))
            sys.exit(1)
    if selected_species not in subfolders:
        sys.stderr.write("A directory for clustering results of the selected species ({}) was not found\n".format(selected_species))
        sys.exit(1)
    r_command = 'Rscript -e \'library(methods); shiny::shinyOptions(in_folder="{}", selected_species="{}"); shiny::runApp("{}", launch.browser = TRUE)\''.format(in_folder, selected_species, current_script_dir)
    os.system(r_command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_folder", type=str, help="Path to the folder with the output of gda_clustering.py")
    parser.add_argument("--selected_species", type=str, help="ID of the selected species, if the gda_out folder contains the results of clustering of more than one species (this has to match a subfolder name in the gda_out folder)", default="")
    args = parser.parse_args()
    main(args.in_folder, args.selected_species)

