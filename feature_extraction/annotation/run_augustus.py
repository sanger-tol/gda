#!/usr/bin/env python3
"""
Script for running Augustus for genome annotation as multiple parallel jobs
This script uses snippets of code adapted from https://github.com/stephenrdoyle/generic_scripts/blob/master/random_workflows/run_augustus_split_by_contigs.sh
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
import os
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import argparse
import tempfile


def main(fasta_path, augustus_folder, hints_gff_path, augustus_species, augustus_chunk_overlap, augustus_chunk_size, threads):
    dry_run_flag = False

    gpf.run_system_command("mkdir -p " + augustus_folder, dry_run=dry_run_flag)

    fasta_basename = fasta_path.split("/")[-1]
    fasta_basename = fasta_basename.split(".")[0]

    os.chdir(augustus_folder)

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        sys.stderr.write("Working directory: " + tmp_dir + "\n")

        splitfasta_command = "splitMfasta.pl {} --outputpath={}".format(fasta_path, tmp_dir)
        gpf.run_system_command(splitfasta_command, dry_run=dry_run_flag)

        rename_fasta_command = "for f in *.split.*; do NAME=`grep \">\" $f`; mv $f ${NAME#>}.fa; done"
        gpf.run_system_command(rename_fasta_command, dry_run=dry_run_flag)

        summarise_acgt_content_command = "summarizeACGTcontent.pl {} > {}/basesummary.out".format(fasta_path, tmp_dir)
        gpf.run_system_command(summarise_acgt_content_command, dry_run=dry_run_flag)

        hints_gff_filename = None
        if hints_gff_path != "NA":
            if os.path.isfile(hints_gff_path):
                gpf.run_system_command("cp {} {}".format(hints_gff_path, tmp_dir), dry_run=dry_run_flag)
                hints_gff_filename = hints_gff_path.split("/")[-1]
            else:
                sys.stderr.write("The hints file for Augustus ({}) was not found\n".format(hints_gff_path))
                sys.exit(1)
        else:
            hints_gff_filename = "empty_hints_file.gff"
            gpf.run_system_command("touch {}/{}".format(tmp_dir, hints_gff_filename), dry_run=dry_run_flag)


        hints_processing_command = "grep \"bases\" basesummary.out | awk -v PWD=$PWD -v HINTS=" + hints_gff_filename + " '{print PWD\"/\"$3\".fa\",PWD\"/\"HINTS,\"1\",$1}' OFS=\"\\t\" > sequences.list"
        gpf.run_system_command(hints_processing_command, dry_run=dry_run_flag)

        split_out_folder_path = tmp_dir + "/augustus_split_out"
        gpf.run_system_command("mkdir -p " + split_out_folder_path, dry_run=dry_run_flag)

        create_joblist_command = "createAugustusJoblist.pl --sequences=sequences.list --wrap=\"#\" --overlap={} --chunksize={} --outputdir=augustus_split_out --joblist=jobs.lst --jobprefix=augsplit --command \"augustus --strand=both --progress=true --gff3=on --codingseq=on --exonnames=on --species={}\"".format(augustus_chunk_overlap, augustus_chunk_size, augustus_species)
        gpf.run_system_command(create_joblist_command, dry_run=dry_run_flag)

        create_augsplit_file_command = "for i in augsplit*; do cat \"${i}\" | grep -v '#' >> run_augsplit; done"
        gpf.run_system_command(create_augsplit_file_command, dry_run=dry_run_flag)

        run_augsplit_command = "parallel --jobs {} < {}/run_augsplit".format(threads, tmp_dir)
        gpf.run_system_command(run_augsplit_command, dry_run=dry_run_flag)

        merged_gff_filename = "{}_augustus_merge.gff".format(fasta_basename)
        filtered_gff_filename = "{}_augustus_merge_filtered.gff".format(fasta_basename)

        merge_gff_command = 'for f in augustus_split_out/*gff; do cat "$f" >> merged_raw.gff; done; cat merged_raw.gff | join_aug_pred.pl > {}'.format(merged_gff_filename)
        gpf.run_system_command(merge_gff_command, dry_run=dry_run_flag)

        filter_gff_command = "grep \"AUGUSTUS\" {} > temp.gff".format(merged_gff_filename)
        gpf.run_system_command(filter_gff_command, dry_run=dry_run_flag)

        sort_gff_command = "echo \"##gff-version 3\" > {0}; bedtools sort -i temp.gff >> {0}".format(filtered_gff_filename)
        gpf.run_system_command(sort_gff_command, dry_run=dry_run_flag)

        copy_gff_command = "cp {} {}".format(filtered_gff_filename, augustus_folder)
        gpf.run_system_command(copy_gff_command, dry_run=dry_run_flag)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("augustus_folder", type=str, help="Folder for running RepeatModeler")
    parser.add_argument("augustus_species", type=str, help="Augustus species")
    parser.add_argument("--hints_gff_path", type=str, help="Optional: path to GFF file with hints for Augustus", default="NA")
    parser.add_argument("--augustus_chunk_overlap", type=int, help="Augustus chunk overlap (in basepairs, default: 12500)", default=12500)
    parser.add_argument("--augustus_chunk_size", type=int, help="Augustus chunk size (in basepairs, default: 750000)", default=750000)
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    args = parser.parse_args()
    main(args.fasta_path, args.augustus_folder, args.hints_gff_path, args.augustus_species, str(args.augustus_chunk_overlap), str(args.augustus_chunk_size), str(args.threads))

