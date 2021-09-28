#!/usr/bin/env python3
"""
Script conversion of .sam file with mapped reads to sorted and indexed .bam file
"""

import general_purpose_functions as gpf
import argparse

def main(sam_file, threads, fasta_path):
    sam_file_extensionless_path = sam_file.split(".sam")[0]
    bam_file = sam_file_extensionless_path + ".bam"
    sorted_bam_file = sam_file_extensionless_path + "_sorted.bam"
    sam_to_bam_command = None
    if fasta_path == "":
        sam_to_bam_command = "samtools view -Sb -@ {} {} > {}".format(threads, sam_file, bam_file)
    else:
        sam_to_bam_command = "samtools view -T {} -Sb -@ {} {} > {}".format(fasta_path, threads, sam_file, bam_file)
    gpf.run_system_command(sam_to_bam_command)
    sort_bam_command = "samtools sort -@ {} {} -o {}".format(threads, bam_file, sorted_bam_file)
    gpf.run_system_command(sort_bam_command)
    index_bam_command = "samtools index " + sorted_bam_file
    gpf.run_system_command(index_bam_command)
    remove_sam_command = "rm " + sam_file
    gpf.run_system_command(remove_sam_command)
    remove_unsorted_bam_command = "rm " + bam_file
    gpf.run_system_command(remove_unsorted_bam_command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sam_file", type=str, help="Path to input SAM file")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    parser.add_argument("--fasta_path", type=str, help="Optional (needed when the reads have been mapped to a genome larger than 4Gb with minimap2): path the FASTA file to which the reads were mapped", default="")
    args = parser.parse_args()
    main(args.sam_file, args.threads, args.fasta_path)