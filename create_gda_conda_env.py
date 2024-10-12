#!/usr/bin/env python3
"""
Script for creating a conda environment for GDA

Prerequisites:
Linux OS (with these standard Linux command line tools: wget, gunzip, tar, sed)
Python3
git
GDA git repository files
g++ 4.9.1 or later
conda
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
import subprocess
import argparse
from glob import glob
from shutil import copyfile
from shutil import which
import platform


def run_system_command(system_command, expected_exit_codes=[0], verbose=True, dry_run=False):
    """
    Executes a system command and checks its exit code
    """
    if verbose == True:
        sys.stderr.write("Executing command: " + system_command + "\n")
    if dry_run == False:
        exit_code = os.system(system_command)
        if exit_code not in expected_exit_codes:
            sys.stderr.write("Error running command: " + system_command + "\n" + "Exited with exit code " + str(exit_code) + "\n")
            sys.exit(1)


def setup_orthomcl(orthomcl_install_path, mcl_path):
    """
    Downloads OrthoMCL 1.4 and edits its orthomcl_module.pm as required
    """
    run_system_command("mkdir -p {}".format(orthomcl_install_path))
    os.chdir(orthomcl_install_path)
    run_system_command("wget http://www.orthomcl.org/common/downloads/software/unsupported/v1.4/ORTHOMCL_V1.4_mcl-02-063.tar -O {}/omcl.tar".format(orthomcl_install_path))
    run_system_command("tar -xvf {}/omcl.tar".format(orthomcl_install_path))
    run_system_command("rm {}/omcl.tar".format(orthomcl_install_path))
    run_system_command("rm {}/mcl-02-063.tar.gz".format(orthomcl_install_path))

    mcl_path_escaped = mcl_path.replace("/", "\\/")

    sed_strings = ['\'s/our .PATH_TO_ORTHOMCL.*=.*/our $PATH_TO_ORTHOMCL = ".\/";/\'',
        '\'s/our .BLASTALL.*=.*/our $BLASTALL = "\/usr\/bin\/blastall";/\'',
        '\'s/our .FORMATDB.*=.*/our $FORMATDB = "\/usr\/bin\/formatdb";/\'',
        '\'s/\/Users\/fengchen\/mcl-02-063\/shmcl\/mcl/{}/g\''.format(mcl_path_escaped)]

    for sed_string in sed_strings:
        sed_command = "sed -i {} {}/ORTHOMCLV1.4/orthomcl_module.pm".format(sed_string, orthomcl_install_path)
        run_system_command(sed_command)

    run_system_command("wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.7.7.tar.gz -O {}/BioPerl-1.7.7.tar.gz".format(orthomcl_install_path))
    run_system_command("tar -xzvf {}/BioPerl-1.7.7.tar.gz".format(orthomcl_install_path))
    run_system_command("rm {}/BioPerl-1.7.7.tar.gz".format(orthomcl_install_path))


def setup_meshclust2(gda_downloads_folder):
    """
    Downloads and compiles MeShClust2 (which also contains a Red binary file)
    """
    os_command = "cd {}; git clone https://github.com/BioinformaticsToolsmith/MeShClust2.git; cd MeShClust2; cd src; mkdir bin; cmake ..; make".format(gda_downloads_folder)
    run_system_command(os_command)


def get_dfam_hmms(dfam_user_folder, gda_downloads_folder, dfam_file_names):
    """
    Triggers the downloading of Dfam database or the copying of it from an user-defined path
    """
    dfam_hmms_folder = gda_downloads_folder + "/dfam_hmm"
    run_system_command("mkdir -p {}".format(dfam_hmms_folder))
    os.chdir(dfam_hmms_folder)
    if dfam_user_folder == "":
        dfam_hmms_base_url = "https://www.dfam.org/releases/Dfam_3.3/infrastructure/dfamscan"
        for dfam_file in dfam_file_names:
            os_command = "wget {}/{}".format(dfam_hmms_base_url, dfam_file)
            run_system_command(os_command)

    else:
        for dfam_file in dfam_file_names:
            src_path = dfam_user_folder + "/" + dfam_file
            dst_path = dfam_hmms_folder + "/" + dfam_file
            copyfile(src_path, dst_path)


def check_user_provided_dfam_folder(dfam_user_folder, dfam_file_names):
    """
    Checks the user-provided Dfam database folder to see if it contains the expected files
    """
    dfam_user_files = os.listdir(dfam_user_folder)
    found_files = [n for n in dfam_file_names if n in dfam_user_files]
    if len(found_files) != len(dfam_file_names):
        not_found_files = [n for n in dfam_file_names if n not in dfam_user_files]
        sys.stderr.write("The expected Dfam database files ({}) were not found in the user-specified folder ({})\n".format(", ".join(not_found_files), dfam_user_folder))
        sys.exit(1)


def check_gda_downloads_folder(gda_downloads_folder):
    """
    Creates the GDA downloads folder if it does not exist yet. Exits with an error message if the folder exists and is not empty
    """
    if os.path.isdir(gda_downloads_folder) == True:
        gda_downloads_folder_files = os.listdir(gda_downloads_folder)
        if len(gda_downloads_folder_files) > 0:
            sys.stderr.write("The folder specified as gda_downloads_folder ({}) already exists and is not empty\n".format(gda_downloads_folder))
            sys.exit(1)
    else:
        run_system_command("mkdir -p " + gda_downloads_folder)


def check_gcc_version(query_version_string, min_gcc_version):
    """
    Parses the string that contains system g++ version and compares it with a float that contains minimum g++ version number. Exits if g++ is too old
    """
    query_version_string = query_version_string.strip()
    allowed_chars = "1234567890."
    gcc_version_reformatted = ""
    dot_encountered_flag = False
    for char in query_version_string:
        if char not in allowed_chars:
            sys.stderr.write("Unrecognised g++ version ({})\n".format(query_version_string))
            sys.exit(1)
        if char != ".":
            gcc_version_reformatted = gcc_version_reformatted + char
        else:
            if dot_encountered_flag == False:
                gcc_version_reformatted = gcc_version_reformatted + "."
                dot_encountered_flag = True
    gcc_version_reformatted = float(gcc_version_reformatted)
    if gcc_version_reformatted < min_gcc_version:
        sys.stderr.write("g++ version {} is too old (should be >= 4.9.1)\n".format(query_version_string))
        sys.exit(1)


def check_prerequisite_software():
    """
    Checks if prerequisite software is installed
    """
    platform_name = platform.system()
    if platform_name != "Linux":
        sys.stderr.write("Warning: this installation script is meant for Linux but your operating system was detected as {}\n".format(platform_name))
    prerequisites_list = ["wget", "gunzip", "tar", "sed", "python3", "git", "conda", "g++"]
    for tool in prerequisites_list:
        if which(tool) == None:
            sys.stderr.write("Prerequisite tool {} was not found\n".format(tool))
            sys.exit(1)
    system_gplusplus_version = subprocess.check_output(["g++", "-dumpfullversion"]).decode("utf-8")
    check_gcc_version(system_gplusplus_version, 4.91)


def main(env_name, gda_downloads_folder, gda_python_scripts_folder, download_dfam_db, dfam_user_folder):

    check_prerequisite_software()

    gda_downloads_folder = gda_downloads_folder.rstrip("/")
    gda_downloads_folder = os.path.abspath(gda_downloads_folder)
    check_gda_downloads_folder(gda_downloads_folder)

    orthomcl_install_folder = gda_downloads_folder + "/orthomcl"
    gda_python_scripts_folder = gda_python_scripts_folder.rstrip("/")
    gda_python_scripts_folder = os.path.abspath(gda_python_scripts_folder)

    meshclust2_bin_folder = gda_downloads_folder + "/MeShClust2/bin"

    if os.path.isdir(gda_python_scripts_folder) == False:
        sys.stderr.write("Directory not found: {}".format(gda_python_scripts_folder))
        sys.exit(1)

    gda_scripts_folder_files = os.listdir(gda_python_scripts_folder)
    if "gda" not in gda_scripts_folder_files or "feature_extraction" not in gda_scripts_folder_files:
        sys.stderr.write("Based on folder contents, the folder {} do not appear to be the base of GDA Python scripts folder\n".format(gda_python_scripts_folder))
        sys.exit(1)

    dfam_file_names = ("Dfam.hmm", "Dfam.hmm.h3f", "Dfam.hmm.h3i", "Dfam.hmm.h3m", "Dfam.hmm.h3p")
    if dfam_user_folder != "":
        check_user_provided_dfam_folder(dfam_user_folder, dfam_file_names)

    packages_to_install = "parallel=20201122 repeatmodeler=2.0.1 diamond=2.0.4 perl-bioperl=1.7.2 pandas=1.1.3 emboss=6.6.0 hisat2=2.2.1 r-base=3.6.3 r-gplots=3.1.1 r-ggplot2=3.3.3 trnascan-se=2.0.6 blast=2.10.1 augustus=3.3.3 minimap2=2.17 samtools=1.10 wgsim=1.0 barrnap=0.9 gffread=0.12.1 mcl=14.137 nextflow=0.30.1 trf=4.09.1 openjdk=8.0.412 mafft=7.475 hmmer=3.3.2 umap-learn=0.4.6 hdbscan=0.8.27 matplotlib=3.3.4 json=0.1.1 numba=0.51.2 spectra=0.0.11 statsmodels=0.12.2 scipy=1.6.2 biopython=1.78 kcounter=0.1.0 liftoff=1.6.1 bedtools=2.30.0" # numba 0.53.0 causes umap to produce errors
    # Packages for the Shiny app: left out for now because they conflict with the rest of the packages
    # r-shiny=1.1.0 r-rjson=0.2.20 r-reshape2=1.4.3 r-gridextra=2.3 r-scales=1.0.0 r-svglite=1.2.3

    channels = ["conda-forge", "bioconda", "jmcmurray"]

    channels_string = "-c " + " -c ".join(channels)

    run_system_command("conda create -y --name {} {} {}".format(env_name, packages_to_install, channels_string))

    conda_base_path = subprocess.check_output(['conda', 'info', '--base'])
    conda_base_path = conda_base_path.decode('utf-8').strip()

    env_path = "{}/envs/{}".format(conda_base_path, env_name)

    mcl_path = env_path + "/bin/mcl"

    gda_feature_extraction_subfolders = glob(gda_python_scripts_folder + "/feature_extraction/*/")
    gda_feature_extraction_subfolders = [n.rstrip("/") for n in gda_feature_extraction_subfolders]

    activate_dir_path = env_path + "/etc/conda/activate.d"
    deactivate_dir_path = env_path + "/etc/conda/deactivate.d"

    run_system_command("mkdir -p {}".format(activate_dir_path))
    run_system_command("mkdir -p {}".format(deactivate_dir_path))

    activate_sh_path = activate_dir_path + "/env_vars.sh"

    with open(activate_sh_path, "w") as f:
        f.write("#!/bin/sh\n")
        f.write("export PATH={}:$PATH\n".format(gda_python_scripts_folder))
        f.write("export PATH={}/feature_extraction:$PATH\n".format(gda_python_scripts_folder))
        f.write("export PATH={}/feature_extraction/cpp/bin:$PATH\n".format(gda_python_scripts_folder))
        for feature_extraction_subfolder in gda_feature_extraction_subfolders:
            f.write("export PATH={}:$PATH\n".format(feature_extraction_subfolder))

        f.write("export PATH={}:$PATH\n".format(meshclust2_bin_folder))

        f.write("export PATH={}/ORTHOMCLV1.4:$PATH\n".format(orthomcl_install_folder))
        f.write("export PERL5LIB={}/ORTHOMCLV1.4:$PERL5LIB\n".format(orthomcl_install_folder))
        f.write("export PERL5LIB={}/BioPerl-1.7.7/lib:$PERL5LIB\n".format(orthomcl_install_folder))
        f.write("export GDA_DOWNLOADS_FOLDER={}\n".format(gda_downloads_folder))
    run_system_command("chmod +x {}".format(activate_sh_path))

    setup_orthomcl(orthomcl_install_folder, mcl_path)
    setup_meshclust2(gda_downloads_folder)
    if download_dfam_db == True:
        get_dfam_hmms(dfam_user_folder, gda_downloads_folder, dfam_file_names)
    sys.stderr.write("Finished setting up the conda environment. It should now be possible to activate the environment with the command 'conda activate {}'\n".format(env_name))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("env_name", type=str, help="Name for the new conda environment (e.g. gda)")
    parser.add_argument("gda_downloads_folder", type=str, help="Path for a new folder where GDA will download some required files that will be installed outside conda")
    parser.add_argument("gda_python_scripts_folder", type=str, help="Path to folder with GDA Python scripts (the scripts should already be present there). This should be the root of the GDA git repository, containing files like 'gda' and 'Singularity'")
    parser.add_argument("--download_dfam_db", dest="download_dfam_db", help="Boolean argument that determines whether Dfam repeat family database is downloaded (the use of Dfam database to classify repeat families is an optional step of the Red+MeShClust2 repeat family detection pipeline)", action="store_true")
    parser.add_argument("--dfam_user_folder", type=str, help="Path to a folder with the Dfam database files if it has already been from Dfam website (https://www.dfam.org/releases/Dfam_3.3/infrastructure/dfamscan) and is locally stored. If a path is not provided but the 'download_dfam_db' flag is set to True, this script will download the database files from the Dfam website", default="")
    args = parser.parse_args()
    main(args.env_name, args.gda_downloads_folder, args.gda_python_scripts_folder, args.download_dfam_db, args.dfam_user_folder)




