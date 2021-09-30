#!/usr/bin/env python3
"""
Script conversion of .sam file with mapped reads to sorted and indexed .bam file
Argument1: path to .sam file
Argument2: number of threads
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

import sys
import general_purpose_functions as gpf

sam_file = sys.argv[1]
threads = sys.argv[2]

sam_file_extensionless_path = sam_file.split(".sam")[0]
bam_file = sam_file_extensionless_path + ".bam"
sorted_bam_file = sam_file_extensionless_path + "_sorted.bam"
sam_to_bam_command = "samtools view -Sb -@ {} {} > {}".format(threads, sam_file, bam_file)
gpf.run_system_command(sam_to_bam_command)
sort_bam_command = "samtools sort -@ {} {} -o {}".format(threads, bam_file, sorted_bam_file)
gpf.run_system_command(sort_bam_command)
index_bam_command = "samtools index " + sorted_bam_file
gpf.run_system_command(index_bam_command)
remove_sam_command = "rm " + sam_file
gpf.run_system_command(remove_sam_command)
remove_unsorted_bam_command = "rm " + bam_file
gpf.run_system_command(remove_unsorted_bam_command)
