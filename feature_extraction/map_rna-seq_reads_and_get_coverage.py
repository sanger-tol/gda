#!/usr/bin/env python3
"""
Script for getting RNA-Seq read coverage for decomposition analysis of genomes
"""

import general_purpose_functions as gpf
import argparse
import os


def main(fasta_path, fastq_1_path, fastq_2_path, mapping_folder, pipeline_output_folder, threads, chunk_size, rf_stranded, library_name, max_intronlen):
    strandedness_flag = ""
    if rf_stranded == True:
        strandedness_flag = "--rna-strandness RF"

    gpf.run_system_command("mkdir -p " + mapping_folder)

    threads = str(threads)
    fasta_filename = fasta_path.split("/")[-1]
    fasta_basename = fasta_filename.split(".")[0]

    if library_name != "":
        library_name = "_" + library_name

    mapping_folder_fasta_path = mapping_folder + "/" + fasta_filename
    output_path_stem = mapping_folder + "/" + fasta_basename + library_name
    sam_file_path = output_path_stem + ".sam"
    sorted_bam_path = output_path_stem + "_sorted.bam"
    summary_file_path = output_path_stem + "_hisat2_summary.txt"
    metrics_file_path = output_path_stem + "_hisat2_metrics.txt"

    samtools_depth_file_path = output_path_stem + "_hisat2_samtools_depth.txt"
    bedgraph_file_path = pipeline_output_folder + "/" +  fasta_basename + library_name + "_hisat2_samtools_depth.bedgraph"

    if os.path.isfile(mapping_folder_fasta_path) == False:
        gpf.run_system_command("cp " + fasta_path + " " + mapping_folder_fasta_path)
        gpf.run_system_command("hisat2-build {} {}".format(mapping_folder_fasta_path, mapping_folder_fasta_path))

    max_intronlen_string = ""
    if max_intronlen != -1:
        max_intronlen_string = "--max-intronlen " + str(max_intronlen) + " "
    hisat_command = "hisat2 -p {} {} -q -x {} -1 {} -2 {} -S {} {} --summary-file {} --met-file {}".format(threads, strandedness_flag, mapping_folder_fasta_path, fastq_1_path, fastq_2_path, sam_file_path, max_intronlen_string, summary_file_path, metrics_file_path)
    gpf.run_system_command(hisat_command)
    sam_to_bam_command = "sam_to_sorted_indexed_bam.py " + sam_file_path + " " + threads
    gpf.run_system_command(sam_to_bam_command)

    samtools_depth_command = "samtools depth -aa {} --reference {} > {}".format(sorted_bam_path, mapping_folder_fasta_path, samtools_depth_file_path)
    gpf.run_system_command(samtools_depth_command)

    depth_to_bedgraph_command = "samtools_depth_to_bedgraph.py {} HISAT2_RNA-Seq_coverage --chunk_size {} > {}".format(samtools_depth_file_path, chunk_size, bedgraph_file_path)
    gpf.run_system_command(depth_to_bedgraph_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("fastq_1_path", type=str, help="Path to read 1 fastq.gz file of RNA-Seq paired end reads")
    parser.add_argument("fastq_2_path", type=str, help="Path to read 2 fastq.gz file of RNA-Seq paired end reads")
    parser.add_argument("mapping_folder", type=str, help="Folder path for RNA-Seq read mapping output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder where to save the output bedgraph file")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--rf_stranded", dest="rf_stranded", help="Boolean argument for strandedness in the mapping of RNA-Seq reads. If this flag is enabled, HISAT2 is run with --rna-strandness RF flag. Otherwise the reads are mapped as unstranded", action="store_true")
    parser.add_argument("--library_name", type=str, help="Optional: RNA-Seq library name (to be used in output file names)", default="")
    parser.add_argument("--max_intronlen", type=int, help="Optional: max_intronlen setting value for HISAT2", default=-1)
    args = parser.parse_args()
    main(args.fasta_path, args.fastq_1_path, args.fastq_2_path, args.mapping_folder, args.pipeline_output_folder, args.threads, args.chunk_size, args.rf_stranded, args.library_name, args.max_intronlen)