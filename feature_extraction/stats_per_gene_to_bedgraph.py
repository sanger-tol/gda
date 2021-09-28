#!/usr/bin/env python3
"""
Script for converting the table of stats per each gene to bedgraph
"""

import general_purpose_functions as gpf
import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict

def coords_overlap(a_start, a_end, b_start, b_end):
    """
    Function for checking if 2 segments overlap
    """
    overlap_flag = not ((a_end < b_start) or (b_end < a_start))
    return overlap_flag


def get_scaff_lengths(fasta_path):
    """
    Extracts sequence lengths from a FASTA file
    """
    fasta_data = gpf.read_fasta_in_chunks(fasta_path)
    seq_len_dict = OrderedDict()
    for header, seq in fasta_data:
        header = header.split()[0]
        seq_len_dict[header] = len(seq)
    return seq_len_dict


def process_scaff(scaff, df, seq_len_dict, chunk_size, feature_types_dict, include_gene_strand_bias):
    """
    Finds the averages in all genomic windows of one scaffold
    """
    scaff_len = seq_len_dict[scaff]
    scaff_df = df[df["scaff"] == scaff]

    segment_start_pos_range = range(0, scaff_len, chunk_size)
    scaff_segment_means_collection = [process_segment(scaff_df, scaff, n, chunk_size, scaff_len, feature_types_dict, include_gene_strand_bias) for n in segment_start_pos_range]
    return scaff_segment_means_collection


def flatten_list_of_lists(l):
    """
    https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    """
    flat_list = [item for sublist in l for item in sublist]
    return flat_list


def output_df_to_bedgraph(output_df, assembly_title, out_folder, feature_types_dict, reference_dataset_title):
    """
    Exports output_df contents in bedgraph format
    """
    input_file_feature_type = feature_types_dict["input_file_feature_type"]
    if reference_dataset_title != "":
        reference_dataset_title += "_"
    for selected_feature in feature_types_dict[input_file_feature_type]:
        out_path = out_folder + "/" + assembly_title + "_" + reference_dataset_title + selected_feature + ".bedgraph"
        with open(out_path, "w") as f:
            feature_title = selected_feature
            bedgraph_header = 'track type=bedGraph name="' + reference_dataset_title + feature_title + '" description="' + reference_dataset_title + feature_title + '" visibility=full color=0,0,255 altColor=0,100,200 priority=20'
            f.write(bedgraph_header + "\n")
            for row_tuple in output_df.iterrows():
                row = row_tuple[1]
                selected_feature_value = str(row[selected_feature])
                if selected_feature_value != "nan":
                    out_line = " ".join([str(row["scaff"]), str(row["segment_start_pos"]), str(row["segment_end_pos"]), selected_feature_value])
                    f.write(out_line + "\n") 


def process_segment(scaff_df, scaff, segment_start_pos, chunk_size, scaff_len, feature_types_dict, include_gene_strand_bias):
    """
    Finds the averages in one genomic window and returns them as a dictionary
    """
    segment_means_dict = None
    segment_end_pos = segment_start_pos + chunk_size
    if segment_end_pos > scaff_len:
        segment_end_pos = scaff_len

    segment_means_dict = {"scaff": scaff, "segment_start_pos": segment_start_pos, "segment_end_pos": segment_end_pos}
    input_file_feature_type = feature_types_dict["input_file_feature_type"]
    for feature in feature_types_dict[input_file_feature_type]:
        segment_means_dict[feature] = None

    if scaff_df.empty == False:
        coords_overlap_vect = np.vectorize(coords_overlap)(segment_start_pos, segment_end_pos, scaff_df["start"], scaff_df["end"])
        scaff_df2 = scaff_df[coords_overlap_vect]
        if scaff_df2.empty == False:
            if input_file_feature_type == "gene_properties":
                if include_gene_strand_bias == True:
                    segment_means_dict["gene_dna_strand_bias"] = abs(scaff_df2["strand"].mean())
                segment_means_dict["gene_length"] = scaff_df2["gene_length"].mean()
                segment_means_dict["exon_count"] = scaff_df2["exon_count"].mean()
                segment_means_dict["gene_average_exon_length"] = scaff_df2["gene_average_exon_length"].mean()
                segment_means_dict["gene_average_intron_length"] = scaff_df2["gene_average_intron_length"].mean()
            elif input_file_feature_type == "OrthoMCL":
                segment_means_dict["paralog_count"] = scaff_df2["paralog_count"].mean()
                segment_means_dict["ortholog_count"] = scaff_df2["ortholog_count"].mean()
                segment_means_dict["protein_conservation_ratio"] = scaff_df2["protein_conservation_ratio"].mean()
                zero_orthologs_count = (scaff_df2["ortholog_count"] == 0).sum()
                segment_means_dict["species_specific_proteins_ratio"] = zero_orthologs_count / scaff_df2.shape[0]
    return segment_means_dict


def main(fasta_path, input_csv_path, input_file_feature_type, out_folder, assembly_title, reference_dataset_title, chunk_size, include_gene_strand_bias):
    df = pd.read_csv(input_csv_path)
    seq_len_dict = get_scaff_lengths(fasta_path)
    gene_features_list = ["gene_length", "exon_count", "gene_average_exon_length", "gene_average_intron_length"]
    if include_gene_strand_bias == True:
        gene_features_list.append("gene_dna_strand_bias")
    feature_types_dict = {"input_file_feature_type": input_file_feature_type, "gene_properties": gene_features_list, "OrthoMCL": ("paralog_count", "ortholog_count", "protein_conservation_ratio", "species_specific_proteins_ratio")}

    segment_means_collection_nested = [process_scaff(n, df, seq_len_dict, chunk_size, feature_types_dict, include_gene_strand_bias) for n in seq_len_dict.keys()]
    segment_means_collection = flatten_list_of_lists(segment_means_collection_nested)

    output_df = pd.DataFrame(segment_means_collection)
    output_df_to_bedgraph(output_df, assembly_title, out_folder, feature_types_dict, reference_dataset_title)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("input_csv_path", type=str, help="Path to input CSV file with feature values per gene")
    parser.add_argument("input_file_feature_type", type=str, help="Type of features in input file (OrthoMCL or gene_properties)", choices={"OrthoMCL", "gene_properties"})
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph files")
    parser.add_argument("assembly_title", type=str, help="Assembly title")
    parser.add_argument("--reference_dataset_title", type=str, default="", help="Optional: reference dataset title. Default:blank")
    parser.add_argument("--chunk_size", type=int, help="Genome chunk size (bp) for generating the bedgraph file using sliding window (default: 5000)", default=5000)
    parser.add_argument("--include_gene_strand_bias", dest="include_gene_strand_bias", help="Boolean argument that determines whether a gene strand bias bedgraph track is produced (based on gene annotations GFF file)", action="store_true")
    args = parser.parse_args()
    main(args.fasta_path, args.input_csv_path, args.input_file_feature_type, args.pipeline_output_folder, args.assembly_title, args.reference_dataset_title, args.chunk_size, args.include_gene_strand_bias)
