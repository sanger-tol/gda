#!/usr/bin/env python3
"""
Script for checking if required software for the GDA feature extraction pipeline is in path
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


