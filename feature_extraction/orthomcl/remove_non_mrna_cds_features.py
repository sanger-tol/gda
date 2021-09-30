#!/usr/bin/env python3
"""
Script for processing a GFF3 file to remove CDS features whose parent feature is something other than 'mRNA'
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

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse


def get_non_mrna_cds_list(gff_data):
    """
    Takes the contents of a GFF3 file as input. Produces a list of CDS features that do not have 'mRNA' as their parent feature
    """
    mrna_list = list()
    cds_dict = dict()
    non_mrna_cds_list = list()

    for line in gff_data:
        if line.startswith("#") == False:

            split_line = line.split("\t")
            feature_name = split_line[2]
            if feature_name == "CDS":
                feature_description = split_line[8]
                cds_id = gpf.spl(feature_description, "ID=", ";")
                cds_parent = gpf.spl(feature_description, ";Parent=", ";")
                if cds_id in cds_dict:
                    existing_cds_parent = cds_dict[cds_id]
                    if existing_cds_parent != cds_parent:
                        sys.stderr.write("The same CDS ID ({}) appears to correspond to multiple mRNA IDs ({} and {})\n".format(cds_id, existing_cds_parent, cds_parent))
                        sys.exit(1)
                else:
                    cds_dict[cds_id] = cds_parent
            elif feature_name == "mRNA":
                feature_description = split_line[8]
                mrna_id = gpf.spl(feature_description, "ID=", ";")
                mrna_list.append(mrna_id)

    for cds_name in cds_dict:
        parent_id = cds_dict[cds_name]
        if parent_id not in mrna_list:
            non_mrna_cds_list.append(cds_name)
            sys.stderr.write("CDS {} does not appear to have an mRNA parent\n".format(cds_name))
    return non_mrna_cds_list


def main(gff_path):
    gff_data = gpf.l(gff_path)
    non_mrna_cds_list = get_non_mrna_cds_list(gff_data)
    for line in gff_data:
        print_line_flag = True
        if line.startswith("#") == False:
            split_line = line.split("\t")
            feature_name = split_line[2]
            if feature_name == "CDS":
                feature_description = split_line[8]
                cds_id = gpf.spl(feature_description, "ID=", ";")
                if cds_id in non_mrna_cds_list:
                    print_line_flag = False
        if print_line_flag == True:
            print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff_path", type=str, help="Path to input GFF3 file")
    args = parser.parse_args()
    main(args.gff_path)



