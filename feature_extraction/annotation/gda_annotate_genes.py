#!/usr/bin/env python3
"""
Master script for running gene annotation scripts for GDA
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

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse
import os


def main(target_assembly_path, pipeline_run_folder, pipeline_output_folder, augustus_species, target_species_id, augustus_chunk_overlap, augustus_chunk_size, reference_assembly_path, reference_gff_path, barrnap_kingdom, threads, chunk_size):
    dry_run_flag = False

    fasta_basename_with_extension = target_assembly_path.split("/")[-1]
    fasta_basename = fasta_basename_with_extension.split(".")[0]

    target_species_id = target_species_id.replace(" ", "_")
    target_species_id = target_species_id.replace("\t", "_")

    gene_annotation_folder = pipeline_run_folder + "/gene_annotation"
    gpf.run_system_command("mkdir -p " + pipeline_run_folder, dry_run=dry_run_flag)
    gpf.run_system_command("mkdir -p " + pipeline_output_folder, dry_run=dry_run_flag)
    gpf.run_system_command("mkdir -p " + gene_annotation_folder, dry_run=dry_run_flag)
    os.chdir(gene_annotation_folder)

    #=== Liftoff ===
    liftoff_folder = gene_annotation_folder + "/liftoff"
    gpf.run_system_command("mkdir -p " + liftoff_folder, dry_run=dry_run_flag)
    os.chdir(liftoff_folder)

    liftoff_output_gff_path = "NA"
    liftoff_augustus_hints_path = "NA"
    if reference_assembly_path != "NA" and reference_gff_path != "NA":
        if os.path.isfile(reference_assembly_path) and os.path.isfile(reference_gff_path):
            liftoff_output_gff_path = liftoff_folder + "/" + fasta_basename + "_liftoff.gff3"
            liftoff_augustus_hints_path = liftoff_folder + "/" + fasta_basename + "_liftoff_augustus_hints.gtf"
            liftoff_stdout_path = liftoff_folder + "/" + fasta_basename + "_liftoff_stdout.txt"

            liftoff_command = "run_liftoff.py {} {} {} {} {} --threads {} > {}".format(target_assembly_path, reference_assembly_path, reference_gff_path, liftoff_output_gff_path, liftoff_augustus_hints_path, threads, liftoff_stdout_path)
            gpf.run_system_command(liftoff_command, dry_run=dry_run_flag)
        else:
            sys.stderr.write("Error: reference file for Liftoff not found:\n")
            if os.path.isfile(reference_assembly_path) == False:
                sys.stderr.write("{}\n".format(reference_assembly_path))
            if os.path.isfile(reference_gff_path) == False:
                sys.stderr.write("{}\n".format(reference_gff_path))
    os.chdir(gene_annotation_folder)

    #=== Augustus ===
    augustus_gff_path = "{}/{}_augustus_merge_filtered.gff".format(gene_annotation_folder, fasta_basename)
    augustus_command = "run_augustus.py {} {} {} --hints_gff_path {} --augustus_chunk_overlap {} --augustus_chunk_size {} --threads {}".format(target_assembly_path, gene_annotation_folder, augustus_species, liftoff_augustus_hints_path, augustus_chunk_overlap, augustus_chunk_size, threads)
    gpf.run_system_command(augustus_command, dry_run=dry_run_flag)

    #=== tRNAscan-SE ===
    trnascan_command = "run_trnascan.py {} {} --threads {}".format(pipeline_run_folder, target_assembly_path, threads)
    gpf.run_system_command(trnascan_command, dry_run=dry_run_flag)
    trnascan_gff_path = gene_annotation_folder + "/trnascan/" + fasta_basename + "_tRNAscan_tRNAs.gff3"

    #=== Barrnap ===
    barrnap_command = "run_barrnap.py {} {} --kingdom {} --threads {}".format(pipeline_run_folder, target_assembly_path, barrnap_kingdom, threads)
    gpf.run_system_command(barrnap_command, dry_run=dry_run_flag)
    barrnap_gff_path = gene_annotation_folder + "/barrnap/" + fasta_basename + "_Barrnap_rRNAs.gff3"

    #=== Combining annotation GFF3 files ===
    combine_annotation_gff_files_command = "combine_annotation_gff_files.py {} {} {} {} {} {}".format(augustus_gff_path, barrnap_gff_path, trnascan_gff_path, target_assembly_path, gene_annotation_folder, target_species_id)
    gpf.run_system_command(combine_annotation_gff_files_command, dry_run=dry_run_flag)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_assembly_path", type=str, help="Path to the assembly that will be annotated")

    parser.add_argument("pipeline_run_folder", type=str, help="Folder path for pipeline output files")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder path for main output files of the pipeline (where bedgraph files will be saved)")
    parser.add_argument("augustus_species", type=str, help="Augustus species")
    parser.add_argument("target_species_id", type=str, help="Species ID of the target species (will be used in feature names in the output GFF file)")

    parser.add_argument("--augustus_chunk_overlap", type=int, help="Augustus chunk overlap (in basepairs, default: 12500)", default=12500)
    parser.add_argument("--augustus_chunk_size", type=int, help="Augustus chunk size (in basepairs, default: 750000)", default=750000)

    parser.add_argument("--reference_assembly_path", type=str, help="Path to the reference assembly for Liftoff, default=NA", default="NA")
    parser.add_argument("--reference_gff_path", type=str, help="Path to the reference GFF file for Liftoff, default=NA", default="NA")
    parser.add_argument("--barrnap_kingdom", type=str, help="Super kingdom of the species for Barrnap. Options: 'bac' (bacteria), 'arc' (archaea), 'euk' (eukaryota), 'mito': metazoan mitochondria. Default: 'euk'", choices=["bac", "arc", "euk", "mito"], default="euk")
    parser.add_argument("--threads", type=int, help="Number of CPU threads (default: 1)", default=1)
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)

    args = parser.parse_args()
    main(args.target_assembly_path, args.pipeline_run_folder, args.pipeline_output_folder, args.augustus_species, args.target_species_id, args.augustus_chunk_overlap, args.augustus_chunk_size, args.reference_assembly_path, args.reference_gff_path, args.barrnap_kingdom, str(args.threads), args.chunk_size)
