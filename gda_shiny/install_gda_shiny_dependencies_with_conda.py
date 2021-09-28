#!/usr/bin/env python3
"""
Script for installing GDA Shiny app's dependencies using conda. The script assumes that the conda environment where you want to install the dependencies has already been created and activated 
Example commands for creating and activating a conda environment for this purpose:
conda create -n gda_env_local r-essentials r-base
conda activate gda_env_local
"""

import os
from shutil import which
import sys
import platform

def run_command(system_command):
    """
    Function for running an OS command
    """
    sys.stderr.write("Executing command: " + system_command + "\n")
    exit_code = os.system(system_command)
    if exit_code != 0:
        sys.stderr.write("Error running command: " + system_command + "\n" + "Exited with exit code " + str(exit_code) + "\n")
        sys.exit(1)


def main():
    if which("conda") == None:
        sys.stderr.write("Conda was not found in path\n")
        sys.exit(1)

    kernel_name = platform.system()
    if kernel_name == "Linux" or kernel_name == "Darwin":
        run_command("conda update -n base conda --yes")
        run_command("conda install -y -c r -c conda-forge r-shiny=1.5.0 r-ggplot2=3.2.1 r-gplots=3.0.3 r-rjson=0.2.20 r-reshape2=1.4.3 r-gridextra=2.3 r-scales=1.0.0 r-svglite=1.2.3")
    else:
        sys.stderr.write("This installation script is meant for Linux and MacOS but your system appears to be {}\n".format(kernel_name))
        sys.exit(1)


if __name__ == "__main__":
    main()
