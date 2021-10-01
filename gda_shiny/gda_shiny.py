#!/usr/bin/env python3
"""
Wrapper script for starting GDA Shiny app
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

