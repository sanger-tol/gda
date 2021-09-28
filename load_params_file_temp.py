#!/usr/bin/env python3
"""
Wrapper script for running GDA
"""
import general_purpose_functions as gpf

def load_param_descriptions(nextflow_param_descriptions_path):
    """
    Loads Nextflow parameter descriptions from a file into a dictionary
    """
    params_dict = dict()
    params_data = gpf.l(nextflow_param_descriptions_path)
    for line in params_data[1:len(params_data)]:
        split_line = line.split(",")
        print(split_line)
        param_name = split_line[0]
        param_type = split_line[1]
        param_description = split_line[2]
        params_dict[param_name] = {"type": param_type, "description": param_description}
    return params_dict

params_dict = load_param_descriptions("/home/ea/sanger_data/tools/python_scripts/genome_decomposition/feature_extraction/gda_nextflow_param_descriptions.dat")
print("params_dict", params_dict)