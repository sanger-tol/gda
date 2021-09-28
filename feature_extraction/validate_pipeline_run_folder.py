#!/usr/bin/env python3
"""
Script for validating the GDA pipeline run folder before running the pipeline
"""

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
