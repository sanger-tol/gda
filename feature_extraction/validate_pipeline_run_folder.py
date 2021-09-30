#!/usr/bin/env python3
"""
Script for validating the GDA pipeline run folder before running the pipeline
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

import os
import sys
import argparse

def main(pipeline_run_folder_path):
    """
    Function for checking if the pipeline run folder is empty before the run starts
    """

    if os.path.isdir(pipeline_run_folder_path) == True:
        pipeline_run_folder_files = os.listdir(pipeline_run_folder_path)
        pipeline_run_folder_files = [n for n in pipeline_run_folder_files if n.endswith(".config") == False]
        pipeline_run_folder_files = [n for n in pipeline_run_folder_files if n not in ("nextflow_trace.txt", ".nextflow.log", ".nextflow", "work")]
        if len(pipeline_run_folder_files) > 0:
            sys.stderr.write("The folder for running the pipeline ({}) already contains files other than Nextflow files ({})\n".format(pipeline_run_folder_path, ", ".join(pipeline_run_folder_files)))
            sys.exit(1)
    else:
        sys.stderr.write("The folder for running the pipeline ({}) does not exist\n".format(pipeline_run_folder_path))
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pipeline_run_folder_path", type=str, help="Path to the empty folder that should serve as the GDA pipeline run folder")
    args = parser.parse_args()
    main(args.pipeline_run_folder_path)
