#!/usr/bin/env python3
"""
Script for concatenating TSV tables (that have been made from bedgraph files) from multiple GDA runs, in order to cluster the data from multiple genomes together. 
Input: paths to TSV files. Output (STDOUT): a concatenated TSV file.
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

import pandas as pd
import os.path
import sys
from collections import OrderedDict
import argparse


def validate_tsv(tsv_file):
    """
    Checks if the input file is a TSV file and has the mandatory columns
    """
    REQUIRED_COLUMNS = ("window", "start", "end", "species", "chromosome")
    if os.path.isfile(tsv_file) == False:
        sys.stderr.write("Input file {} was not found\n".format(tsv_file))
        sys.exit(1)
    try:
        df = pd.read_csv(tsv_file, sep="\t")
    except pd.errors.ParserError as e:
        sys.stderr.write("Failed to parse the input file {} as a TSV (tab separated values) file\n".format(tsv_file))
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)
    except ValueError as e:
        sys.stderr.write("Failed to parse the input file {} as a TSV (tab separated values) file\n".format(tsv_file))
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)
    except:
        sys.stderr.write("Failed to parse the input file {} as a TSV (tab separated values) file\n".format(tsv_file))
        sys.stderr.write("Unexpected error: " + str(sys.exc_info()[0]) + "\n")
        sys.exit(1)
    df_colnames = df.columns
    
    for required_col in REQUIRED_COLUMNS:
        if required_col not in df_colnames:
            sys.stderr.write("Required column with the name '{}' was not found in the input TSV file {}\n".format(required_col, tsv_file))
            sys.exit(1)


def detect_window_size(df, tsv_file):
    """
    Detects window size in the input TSV file
    """
    start_coords = df["start"].tolist()
    diff_list = list()
    previous_start_coord = 0
    for start_coord in start_coords:
        if start_coord > previous_start_coord:
            diff = start_coord - previous_start_coord
            diff_list.append(diff)
        previous_start_coord = start_coord
    diff_set = set(diff_list)
    window_size = None
    if len(diff_set) > 1:
        sys.stderr.write("More than one window step size ({}) was found in the file {}\n".format(",".join(diff_set), tsv_file))
        sys.exit(1)
    else:
        window_size = list(diff_set)[0]
    return window_size


def check_for_noncommon_columns(df_dict):
    """
    Finds variable columns that are not shared between the input TSV files.
    These columns will be left out from the output TSV file
    """
    common_cols = list(set.intersection(*(set(df.columns) for df in list(df_dict.values()))))
    for tsv_file in df_dict:
        noncommon_cols = [n for n in df_dict[tsv_file].columns if n not in common_cols]
        if len(noncommon_cols) > 0:
            sys.stderr.write("Warning: the following column(s) from the file {} will not be used, as they do not occur in the other TSV file(s): {}\n".format(tsv_file, ", ".join(noncommon_cols)))


def check_for_duplicate_entries(df_dict):
    """
    The values in the columns 'species', 'chromosome' and 'window' have to be different in each input TSV file.
    This function checks if they are different and exits the script if two or more TSV files share the same value
    """
    col_labels = ["species", "chromosome", "window"]
    for col_label in col_labels:
        entries_collection = list()
        for tsv_file in df_dict:
            entries_set = set(df_dict[tsv_file][col_label].tolist())
            for col_value in entries_set:
                if col_value in entries_collection:
                    sys.stderr.write("Value {} appears in the '{}' column in more than one input file\n".format(col_value, col_label))
                    sys.exit(1)
                else:
                    entries_collection.append(col_value)


def check_window_sizes(df_dict):
    """
    The input TSV files need to have the same window size.
    This function checks if the window sizes are equal and exits if they are not
    """
    window_sizes = dict()
    for tsv_file in df_dict:
        window_size = detect_window_size(df_dict[tsv_file], tsv_file)
        window_sizes[tsv_file] = window_size
    if len(set(window_sizes.values())) > 1:
        sys.stderr.write("The input TSV files appear to have different window sizes.\n")
        for tsv_file in window_sizes:
            sys.stderr.write("{}: window size {} bp\n".format(tsv_file, window_sizes[tsv_file]))
        sys.exit(1)


def main(tsv_files):

    for tsv_file in tsv_files:
        validate_tsv(tsv_file)

    df_dict = OrderedDict()
    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep="\t")
        if tsv_file not in df_dict:
            df_dict[tsv_file] = df
        else:
            sys.stderr.write("The input file {} appears in the input more than once\n".format(tsv_file))
            sys.exit(1)

    check_for_noncommon_columns(df_dict)
    check_window_sizes(df_dict)
    check_for_duplicate_entries(df_dict)
    
    concat_df = pd.concat(list(df_dict.values()), join="inner", ignore_index=True, sort=False)
    concat_df.to_csv(sys.stdout, index=False, sep="\t", header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv_files", type=str, nargs="+", default=[], help="Path(s) to the TSV file(s)")
    args = parser.parse_args()
    main(args.tsv_files)