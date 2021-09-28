#!/usr/bin/env python3
"""
Script for extracting mRNA, pseudogene, tRNA and rRNA features from GFF3 file to convert them into bedgraph format
"""

import general_purpose_functions as gpf
import argparse


def main(gff_path, fasta_path, out_folder, chunk_size, mRNA_tag, pseudogene_tag, tRNA_tag, rRNA_tag, custom_gff_tags):

    fasta_basename = fasta_path.split("/")[-1]
    fasta_basename = fasta_basename.split(".")[0]

    gpf.run_system_command("mkdir -p " + out_folder)

    selected_tags = [mRNA_tag, pseudogene_tag, tRNA_tag, rRNA_tag]

    if custom_gff_tags != "" and custom_gff_tags != "NA":
        split_custom_gff_tags = custom_gff_tags.split(",")
        split_custom_gff_tags = [n for n in split_custom_gff_tags if n != ""]
        split_custom_gff_tags = [n.strip() for n in split_custom_gff_tags]
        split_custom_gff_tags = list(set(split_custom_gff_tags))
        selected_tags.extend(split_custom_gff_tags)


    for selected_tag in selected_tags:
        if selected_tag != "NA":
            bedgraph_path = out_folder + "/" + fasta_basename + "_" + selected_tag + ".bedgraph"
            os_command = "gff_features_to_bedgraph.py {} {} {}_annotations --assembly_fasta_path {} --chunk_size {} > {}".format(gff_path, selected_tag, selected_tag, fasta_path, str(chunk_size), bedgraph_path)
            gpf.run_system_command(os_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff_path", type=str, help="Path to input GFF3 file")
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder for outputting bedgraph files")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--mRNA_tag", type=str, help="GFF feature tag for mRNA (default: mRNA)", default="mRNA")
    parser.add_argument("--pseudogene_tag", type=str, help="GFF feature tag for pseudogenes (default: pseudogene)", default="pseudogene")
    parser.add_argument("--tRNA_tag", type=str, help="GFF feature tag for tRNA (default: tRNA)", default="tRNA")
    parser.add_argument("--rRNA_tag", type=str, help="GFF feature tag for rRNA (default: rRNA)", default="rRNA")
    parser.add_argument("--custom_gff_tags", type=str, help="Optional: comma separated string of custom tags for features to extract from the annotations GFF (example: five_prime_UTR,three_prime_UTR,miRNA)", default="")
    args = parser.parse_args()
    main(args.gff_path, args.fasta_path, args.out_folder, args.chunk_size, args.mRNA_tag, args.pseudogene_tag, args.tRNA_tag, args.rRNA_tag, args.custom_gff_tags)
    

