#!/usr/bin/env python3
"""
Script for running Red and MeShClust2 to detect repeat families
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

import general_purpose_functions as gpf
from collections import defaultdict
from collections import OrderedDict
import os
import tempfile
import argparse
import sys


class RepeatCluster:
    """
    Class for storing information on a MeShClust2 repeat cluster
    """
    def __init__(self, id):
        self.id = id
        self.repeats = list()

    def get_repeat_count(self):
        return len(self.repeats)

    def get_average_repeat_length(self):
        """
        Returns average repeat length in the cluster
        """
        lengths_sum = sum([n.length for n in self.repeats])
        average_len = lengths_sum / len(self.repeats)
        return average_len

    def get_repeat_seq(self, repeats_fasta_dict):
        """
        Returns the names and sequences of the repeats in the cluster as an OrderedDict object
        """
        out_dict = OrderedDict()
        for repeat in self.repeats:
            repeat_suffix = "_no_star"
            if repeat.star == True:
                repeat_suffix = "_star"
            out_dict[repeat.name + repeat_suffix] = repeats_fasta_dict[repeat.name]
        return out_dict


class RepeatSeq:
    """
    Class for storing the values related to an individual repeat sequence
    """
    def __init__(self, id, length, star, start_coord, end_coord, scaff):
        self.name = id
        self.length = length
        self.star = star
        self.start_coord = start_coord
        self.end_coord = end_coord
        self.scaff = scaff


def read_repeat_coords_as_dict(red_repeat_coords_file_path):
    """
    Input: path to the repeat coordinates output file of Red
    Output: dictionary where keys are scaffold names and values are a list of strings with repeat coordinates in the corresponding scaffold
    """
    repeat_coords_dict = defaultdict(list)
    repeat_coords_data = gpf.ll(red_repeat_coords_file_path)
    for line in repeat_coords_data:
        line = line.lstrip(">")
        split_line = line.split(":")
        repeat_coords_dict[split_line[0]].append(split_line[1])
    return repeat_coords_dict


def extract_repeat_sequences(assembly_fasta_path, repeat_coords_dict, out_path):
    """
    Extracts the sequences of repeats detected by Red and writes them in a FASTA file
    """
    fasta_data = gpf.read_fasta_in_chunks(assembly_fasta_path)
    with open(out_path, "w") as f:
        for header, seq in fasta_data:
            if header in repeat_coords_dict:
                repeat_coords_list = repeat_coords_dict[header]
                for repeat_coords_string in repeat_coords_list:
                    split_coords = repeat_coords_string.split("-")
                    start_coord = int(split_coords[0])
                    end_coord = int(split_coords[1])
                    repeat_seq = seq[start_coord:end_coord]
                    f.write(">{}_{}\n".format(header, repeat_coords_string))
                    repeat_seq_chunks = gpf.string_to_chunks(repeat_seq, 80)
                    for line in repeat_seq_chunks:
                        f.write(line + "\n")


def run_red(assembly_fasta_path, red_input_folder, red_output_folder, threads):
    """
    Runs Red to detect repeats in an assembly FASTA file
    """
    assembly_fasta_filename = assembly_fasta_path.split("/")[-1]
    renamed_fasta_filename = None
    truncated_filename = None
    if "." in assembly_fasta_filename:
        split_filename = assembly_fasta_filename.split(".")
        truncated_filename = ".".join(split_filename[0:len(split_filename) - 1])
        renamed_fasta_filename = truncated_filename + ".fa"
    else:
        renamed_fasta_filename = assembly_fasta_filename + ".fa"

    renamed_fasta_path = red_input_folder + "/" + renamed_fasta_filename

    gpf.run_system_command("ln -s {} {}".format(assembly_fasta_path, renamed_fasta_path))
    red_command = "Red -gnm {} -rpt {}".format(red_input_folder, red_output_folder)
    gpf.run_system_command(red_command)

    output_file_stem = renamed_fasta_filename[0:len(renamed_fasta_filename) - 3]
    red_repeat_coords_file_path = red_output_folder + "/" + output_file_stem + ".rpt"

    repeat_coords_dict = read_repeat_coords_as_dict(red_repeat_coords_file_path)

    red_repeat_seq_path = red_output_folder + "/red_repeats.fasta"
    extract_repeat_sequences(assembly_fasta_path, repeat_coords_dict, red_repeat_seq_path)


def load_meshclust2_clusters_file(meshclust_outfile_path):
    """
    Loads MeShClust2 output file contents as a dictionary of RepeatCluster objects
    """
    meshclust_data = gpf.l(meshclust_outfile_path)

    clusters_dict = dict()
    current_cluster_name = None
    for line in meshclust_data:
        line = line.strip()
        if line.startswith(">"):
            if current_cluster_name is not None:
                clusters_dict[current_cluster_name] = current_cluster
            current_cluster_name = "cluster_" + line.lstrip(">Cluster ")
            current_cluster = RepeatCluster(current_cluster_name)

        else:
            repeat_len = int(gpf.spl(line, "\t", "nt"))
            repeat_id = gpf.spl(line, ">", "...")
            repeat_star = False
            if line.endswith("*"):
                repeat_star = True
            split_line = line.split("_")
            repeat_start_end = split_line[-1]
            split_repeat_start_end = repeat_start_end.split("-")
            repeat_start_coord = int(split_repeat_start_end[0])
            repeat_end_coord = split_repeat_start_end[1]
            repeat_end_coord = int(repeat_end_coord.split("...")[0])
            repeat_scaff = line.split(">")[1]
            repeat_scaff = repeat_scaff.split("_")
            repeat_scaff = "_".join(repeat_scaff[0:len(repeat_scaff) - 1])

            current_repeat = RepeatSeq(repeat_id, repeat_len, repeat_star, repeat_start_coord, repeat_end_coord, repeat_scaff)
            current_cluster.repeats.append(current_repeat)
    if len(current_cluster.repeats) > 0:
        clusters_dict[current_cluster_name] = current_cluster
    return clusters_dict


def export_meshclust2_results_as_gff(clusters_dict, out_path):
    """
    Exports the locations of repeats in MeShClust2 repeat clusters as a GFF3 file
    """
    f, tempfile_name = tempfile.mkstemp()
    with open(tempfile_name, "w") as f:
        f.write("##gff-version 3\n")
        for cluster_name, cluster in clusters_dict.items():
            for repeat in cluster.repeats:
                gff_line = "{}\tMeShClust2\t{}\t{}\t{}\t.\t+\t.\tID={};repeat_cluster={}".format(repeat.scaff, cluster_name, str(repeat.start_coord + 1), str(repeat.end_coord), repeat.name, cluster_name)
                f.write(gff_line + "\n")
    gff_sort_command = "echo \'##gff-version 3\' > {}; bedtools sort -i {} >> {}".format(out_path, tempfile_name, out_path)
    gpf.run_system_command(gff_sort_command)


def export_repeat_families_as_fasta(red_repeat_seq_path, clusters_dict, min_repeats_per_cluster, repeat_families_outfolder):
    """
    Exports the sequences of repeat clusters detected using Red + MeShClust2 as FASTA files
    """
    repeats_fasta_data = gpf.read_fasta_in_chunks(red_repeat_seq_path)
    repeats_fasta_dict = dict()
    for header, seq in repeats_fasta_data:
        repeats_fasta_dict[header] = seq

    for cluster_name in clusters_dict:
        repeatcluster = clusters_dict[cluster_name]
        cluster_repeat_count = repeatcluster.get_repeat_count()
        if cluster_repeat_count >= min_repeats_per_cluster:
            repeats_dict = repeatcluster.get_repeat_seq(repeats_fasta_dict)
            repeats_file_outpath = repeat_families_outfolder + "/" + cluster_name + ".fa"

            with open(repeats_file_outpath, "w") as f:
                for repeat_name, repeat_seq in repeats_dict.items():
                    f.write(">" + repeat_name + "\n")
                    repeat_seq_chunks = gpf.string_to_chunks(repeat_seq, 80)
                    for line in repeat_seq_chunks:
                        f.write(line + "\n")


def run_meshclust2(red_output_folder, pipeline_output_folder, min_repeats_per_cluster, threads, meshclust2_id, assembly_fasta_path, assembly_name_stem):
    """
    Runs MeShClust2 and extracts the sequences of repeat clusters as FASTA
    """
    os.chdir(red_output_folder)
    repeat_families_outfolder = red_output_folder + "/repeat_families"
    gpf.run_system_command("mkdir -p {}".format(repeat_families_outfolder))

    red_repeat_seq_path = red_output_folder + "/red_repeats.fasta"

    meshclust2_command = "meshclust2 --id {} --threads {} {}".format(str(meshclust2_id), str(threads), red_repeat_seq_path)
    gpf.run_system_command(meshclust2_command)

    meshclust_outfile_path = red_output_folder + "/output.clstr"

    clusters_dict = load_meshclust2_clusters_file(meshclust_outfile_path)
    export_repeat_families_as_fasta(red_repeat_seq_path, clusters_dict, min_repeats_per_cluster, repeat_families_outfolder)

    meshclust2_gff_path = red_output_folder + "/meshclust2_repeat_clusters.gff3"
    export_meshclust2_results_as_gff(clusters_dict, meshclust2_gff_path)

    clusters_bedgraph_folder = pipeline_output_folder + "/meshclust2_complex_repeats"
    gpf.run_system_command("mkdir -p {}".format(clusters_bedgraph_folder))

    gff_to_bedgraph_command = "repeatmasker_gff_to_bedgraph.py {} {} {} complex_repeats --occurrences_cutoff 5 --assembly_fasta_path {}".format(meshclust2_gff_path, clusters_bedgraph_folder, assembly_name_stem, assembly_fasta_path)
    gpf.run_system_command(gff_to_bedgraph_command)

    meshclust2_repeats_sum_path = pipeline_output_folder + "/{}_meshclust2_repeats_sum.bedgraph".format(assembly_name_stem)
    sum_bedgraph_command = "sum_simple_or_complex_repeat_tracks.py {} {} {} complex_repeats > {}".format(clusters_bedgraph_folder, assembly_fasta_path, assembly_name_stem, meshclust2_repeats_sum_path)
    gpf.run_system_command(sum_bedgraph_command)



def get_meshclust2_consensus_seq(repeat_families_folder, alignments_folder, mafft_consensus_folder, threads):
    """
    Runs MAFFT and hmmemit to align sequences in each MeShClust2 cluster and to derive the consensus sequence for each cluster
    """
    gpf.run_system_command("mkdir -p " + alignments_folder)
    gpf.run_system_command("mkdir -p " + mafft_consensus_folder)
    repeat_families_files = [n for n in os.listdir(repeat_families_folder) if n.endswith(".fa")]
    for repeats_file in repeat_families_files:
        file_name_stem = repeats_file.rstrip(".fa")
        alignment_file_name = file_name_stem + "_mafft.fa"
        repeats_file_path = repeat_families_folder + "/" + repeats_file
        mafft_alignment_path = alignments_folder + "/" + alignment_file_name
        mafft_command = "mafft --thread {} --adjustdirection --reorder {} > {}".format(str(threads), repeats_file_path, mafft_alignment_path)
        gpf.run_system_command(mafft_command)
        hmm_file_path = mafft_alignment_path.rstrip("_mafft.fa") + "_mafft.hmm"
        consensus_file_path = mafft_consensus_folder + "/" + file_name_stem + "_mafft_consensus.fa"
        hmmbuild_command = "hmmbuild --dna {} {}".format(hmm_file_path, mafft_alignment_path)
        gpf.run_system_command(hmmbuild_command)
        hmmemit_command = "hmmemit -c {} > {}".format(hmm_file_path, consensus_file_path)
        gpf.run_system_command(hmmemit_command)


def get_dfam_db_path():
    """
    Looks for the path to the Dfam repeats HMM database based on the GDA_DOWNLOADS_FOLDER environmental variable
    Returns None if it is not found
    Unzips the database and runs hmmpress if a gzipped database file is found
    """
    dfam_db_path = None
    if "GDA_DOWNLOADS_FOLDER" in os.environ:
        gda_downloads_folder = os.environ["GDA_DOWNLOADS_FOLDER"]
        if os.path.isdir(gda_downloads_folder) == True:
            dfam_hmms_folder = gda_downloads_folder + "/dfam_hmm"
            if os.path.isdir(dfam_hmms_folder) == True:
                dfam_folder_files = os.listdir(dfam_hmms_folder)
                if "Dfam.hmm" in dfam_folder_files:
                    dfam_db_path = dfam_hmms_folder + "/Dfam.hmm"
                else:
                    sys.stderr.write("GDA Dfam database was not found in the expected folder ({})\n".format(dfam_hmms_folder))
            else:
                sys.stderr.write("GDA Dfam database folder ({}) was not found\n".format(dfam_hmms_folder))
        else:
            sys.stderr.write("GDA downloads folder ({}) was not found\n".format(gda_downloads_folder))
    else:
        sys.stderr.write("GDA_DOWNLOADS_FOLDER environmental variable was not found\n")
    return dfam_db_path



def main(assembly_fasta_path, red_run_folder, pipeline_output_folder, min_repeats_per_cluster, meshclust2_id, chunk_size, threads):

    assembly_name_stem = assembly_fasta_path.split("/")[-1]
    split_assembly_name_stem = assembly_name_stem.split(".")
    assembly_name_stem = ".".join(split_assembly_name_stem[0:len(split_assembly_name_stem) - 1])

    if meshclust2_id > 1 or meshclust2_id < 0.1:
        sys.stderr.write("Invalid value for meshclust2_id ({}). Should be below 1 and above 0.1\n".format(str(meshclust2_id)))
        sys.exit(1)

    gpf.run_system_command("mkdir -p {}".format(red_run_folder))
    red_input_folder = red_run_folder + "/input"
    red_output_folder = red_run_folder + "/output"

    gpf.run_system_command("mkdir -p {}".format(red_input_folder))
    gpf.run_system_command("mkdir -p {}".format(red_output_folder))

    run_red(assembly_fasta_path, red_input_folder, red_output_folder, threads)

    run_meshclust2(red_output_folder, pipeline_output_folder, min_repeats_per_cluster, threads, meshclust2_id, assembly_fasta_path, assembly_name_stem)

    repeat_families_outfolder = red_output_folder + "/repeat_families"
    alignments_folder = red_output_folder + "/repeat_families_mafft"
    gpf.run_system_command("mkdir -p {}".format(alignments_folder))

    mafft_consensus_folder = red_output_folder + "/repeat_families_consensus"
    gpf.run_system_command("mkdir -p {}".format(mafft_consensus_folder))

    get_meshclust2_consensus_seq(repeat_families_outfolder, alignments_folder, mafft_consensus_folder, threads)
    collected_consensus_path = red_output_folder + "/" + assembly_name_stem + "_meshclust2_consensus.fa"
    cat_command = "cat {}/* > {}".format(mafft_consensus_folder, collected_consensus_path)
    gpf.run_system_command(cat_command)

    dfam_db_path = get_dfam_db_path()
    if dfam_db_path is not None:
        dfam_scan_outpath = red_output_folder + "/" + assembly_name_stem + "_meshclust2_clusters_dfam.txt"
        hmmscan_command = "hmmscan --tblout {} --cpu {} {} {}".format(dfam_scan_outpath, str(threads), dfam_db_path, collected_consensus_path)
        gpf.run_system_command(hmmscan_command)
    else:
        sys.stderr.write("Skipping running hmmscan of repeat consensus sequences against the Dfam database, as the database was not found\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly_fasta_path", type=str, help="Path to assembly FASTA file")
    parser.add_argument("red_run_folder", type=str, help="Path to folder where Red will be run")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder path for main output files of the pipeline (where bedgraph files will be saved)")
    parser.add_argument("--min_repeats_per_cluster", type=int, help="Minimum number of repeats in a repeat cluster (default: 5)", default=5)
    parser.add_argument("--meshclust2_id", type=int, help="Identity threshold for repeats when clustering with MeShClust2 (default: 0.95)", default=0.95)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--threads", type=int, help="Number of threads for running BLAST (default: 1)", default=1)
    args = parser.parse_args()
    main(args.assembly_fasta_path, args.red_run_folder, args.pipeline_output_folder, args.min_repeats_per_cluster, args.meshclust2_id, args.chunk_size, args.threads)





