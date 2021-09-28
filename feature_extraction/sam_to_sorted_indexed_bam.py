#!/usr/bin/env python3
"""
Script conversion of .sam file with mapped reads to sorted and indexed .bam file
Argument1: path to .sam file
Argument2: number of threads
"""

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