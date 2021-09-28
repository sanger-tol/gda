#!/usr/bin/env python3
"""
Script for validating the nextflow.config file of GDA
"""

import general_purpose_functions as gpf
from collections import OrderedDict
import sys
import os
import pandas as pd
import numpy as np
from validate_input_files import UnicodeFileChecker
import argparse


class Param:
    """
    Base class for parameters
    """
    def __init__(self, name, value):
        self.name = name
        self.value = value


class IntParam(Param):
    """
    Class for validating integer parameters
    """
    def __init__(self, name, value, min_val, max_val):
        super().__init__(name, value)
        self.min_val = min_val
        self.max_val = max_val

    def validate_param(self):
        try:
            self.value = int(self.value)
        except ValueError:
            sys.stderr.write("Failed to convert the value of the parameter {} ({}) to integer\n")
            sys.exit(1)
        if np.isnan(self.min_val) == False:
            if self.value < self.min_val:
                sys.stderr.write("The selected value of the parameter {} ({}) is below the minimum value for this parameter ({})\n".format(self.name, str(self.value), str(self.min_val)))
                sys.exit(1)
        if np.isnan(self.max_val) == False:
            if self.value > self.max_val:
                sys.stderr.write("The selected value of the parameter {} ({}) is above the maximum allowed value for this parameter ({})\n".format(self.name, str(self.value), str(self.max_val)))
                sys.exit(1)


class MultipleChoiceParam(Param):
    """
    Class for validating multiple choice parameters
    """
    def __init__(self, name, value, choices_string):
        super().__init__(name, value)
        self.choices_string = choices_string

    def validate_param(self):
        choices = self.choices_string.split("/")
        if self.value not in choices:
            choices_output_string = "[" + ",".join(choices) + "]"
            sys.stderr.write("The selected value of the parameter {} ({}) is not among the allowed choices for this variable {}\n".format(self.name, self.value, choices_output_string))
            sys.exit(1)


class BoolParam(Param):
    """
    Class for validating boolean parameters
    """
    def __init__(self, name, value):
        super().__init__(name, value)

    def validate_param(self):
        if self.value != "true" and self.value != "false":
            sys.stderr.write("The selected value of the parameter {} ({}) is not among the allowed choices for this variable (either 'true' or 'false')\n".format(self.name, self.value))
            sys.exit(1)


class StringListParam(Param):
    """
    Class for validating parameters that are a list of strings
    """
    def __init__(self, name, value):
        super().__init__(name, value)

    def validate_param(self):
        if "," in self.value:
            items = self.value.split(",")
            items = [n.strip() for n in items]
            items = [n for n in items if n != ""]
            items2 = list(set(items))
            if len(items2) != len(items):
                sys.stderr.write("Warning: duplicate values found in values of the parameter {} ({})\n".format(self.name, self.value))


class LimitedStringParam(Param):
    """
    Class for validating parameters that are strings with only a limited set of allowed characters
    """
    def __init__(self, name, value):
        super().__init__(name, value)

    def validate_param(self):
        allowed_chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-."
        bad_chars = [n for n in self.value if n not in allowed_chars]
        if len(bad_chars) > 0:
            sys.stderr.write("The selected value of the parameter {} ('{}') contains characters that are not allowed for this variable ('{}'). The allowed characters are these: '{}'\n".format(self.name, self.value, "".join(bad_chars), allowed_chars))
            sys.exit(1)


def remove_trace_field(config_data):
    """
    Input: Nextflow config file loaded as a list of strings
    Output: A list that is like the input list but without the lines that define the Nextflow trace variable, as 
        this part of the config file does not need to be checked in the same way as the rest of the file
    """
    out_list = list()
    trace_line_flag = False
    for line in config_data:
        if line == "trace {":
            trace_line_flag = True
        if trace_line_flag == False:
            out_list.append(line)
        if line == "}":
            trace_line_flag = False
    return out_list
        
        
def main(nextflow_config_path):
    param_descriptions_path = os.path.dirname(os.path.realpath(__file__)) + "/gda_nextflow_param_descriptions.dat"
    gpf.check_if_file_exists(param_descriptions_path)
    df = pd.read_csv(param_descriptions_path)

    gpf.check_if_file_exists(nextflow_config_path)
    config_data = gpf.l(nextflow_config_path)
    unicode_check_file = UnicodeFileChecker(nextflow_config_path)
    unicode_check_file.check_if_file_is_unicode()
    config_data = remove_trace_field(config_data)
    config_data = [n for n in config_data if n.startswith("//") == False and len(n) > 0 and "=" in n]

    config_dict = OrderedDict()
    for line in config_data:
        split_line = line.split("=")
        config_dict[split_line[0].strip()] = split_line[1].strip()
    
    for param_name, param_value in config_dict.items():
        if param_name not in df.param.values:
            sys.stderr.write("Warning: parameter '{}' is not recognised and its value will not be used\n".format(param_name))
        else:
            df_pos = df.loc[df["param"] == param_name].index[0]
            df_row = df.loc[[df_pos]]
            param_name = df_row["param"].values[0]
            param_type = df_row["param_type"].values[0]
            if isinstance(param_value, str):
                param_value = param_value.replace("\"", "")
            if df_row["required"].values[0] == "y":
                if str(param_value) == "NA" or str(param_type).strip() == "":
                    sys.stderr.write("Parameter {} is mandatory with no default value and needs to be assigned a value\n".format(param_name))
                    sys.exit(1)

            if param_type == "int":
                min_val = df_row["min"].values[0]
                max_val = df_row["max"].values[0]
                p = IntParam(param_name, param_value, min_val, max_val)
                p.validate_param()
            elif param_type == "multiple_choice_string" and param_value != "NA":
                choices_string = str(df_row["choices"].values[0])
                p = MultipleChoiceParam(param_name, str(param_value), choices_string)
                p.validate_param()
            elif param_type == "bool":
                p = BoolParam(param_name, str(param_value))
                p.validate_param()
            elif param_type == "stringlist":
                p = StringListParam(param_name, str(param_value))
                p.validate_param()
            elif param_type == "limited_string":
                p = LimitedStringParam(param_name, str(param_value))
                p.validate_param()

                if config_dict["params.run_gene_annotation_pipeline"] == "true" and config_dict["params.augustus_species"] == "NA":
                    sys.stderr.write("Parameter 'params.augustus_species' needs to be assigned a value when 'params.run_gene_annotation_pipeline' has been set to 'true'\n")
                    sys.exit(1)

                if (config_dict["params.reference_assembly_path"] == "NA" and config_dict["params.reference_gff_path"] != "NA") or (config_dict["params.reference_assembly_path"] != "NA" and config_dict["params.reference_gff_path"] == "NA"):
                    sys.stderr.write("To use Liftoff for annotation transfer, 'params.reference_assembly_path' and 'params.reference_gff_path' both need to be assigned values\n")
                    sys.exit(1)

                if (config_dict["params.reference_assembly_path"] == "NA" and config_dict["params.reference_gff_path"] != "NA") or (config_dict["params.reference_assembly_path"] != "NA" and config_dict["params.reference_gff_path"] == "NA"):
                    sys.stderr.write("To use Liftoff for annotation transfer, 'params.reference_assembly_path' and 'params.reference_gff_path' both need to be assigned values\n")
                    sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nextflow_config_path", type=str, help="Path to the nextflow.config file for the Nextflow run")
    args = parser.parse_args()
    main(args.nextflow_config_path)