#!/usr/bin/env python3
"""
Script for checking if required software for the GDA feature extraction pipeline is in path
"""
import platform
from shutil import which
import sys

platform_name = platform.system()
if platform_name != "Linux":
    sys.stderr.write("Warning: this pipeline has been developed and tested under Linux but your operating system was detected as {}\n".format(platform_name))

prerequisites_list = ["augustus",
    "barrnap",
    "blastn",
    "BuildDatabase",
    "conda",
    "diamond",
    "einverted",
    "gffread",
    "hisat2",
    "hmmscan",
    "java",
    "mafft",
    "mcl",
    "meshclust2",
    "minimap2",
    "orthomcl.pl",
    "parallel",
    "python3",
    "R",
    "Red",
    "RepeatMasker",
    "RepeatModeler",
    "samtools",
    "trf",
    "tRNAscan-SE",
    "wgsim"]


for tool in prerequisites_list:
    if which(tool) == None:
        sys.stderr.write("A prerequisite tool for running the GDA genomic feature extraction pipeline ({}) was not found in path\n".format(tool))
        sys.exit(1)


