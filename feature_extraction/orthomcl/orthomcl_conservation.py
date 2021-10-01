#!/usr/bin/env python3
"""
Script for converting OrthoMCL results into a table of paralog counts, ortholog counts and conservation ratio
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
import pandas as pd


def get_gene_ids_dict(annotations_gff_path, query_gff_feature):
    """
    Input: 1) path to the GFF gene annotations file of the species of interest, 2) feature with mRNA gene names (default: 'mRNA')
    Output: dictionary where the keys are the gene IDs. Each entry in the dictionary for a gene ID contains another dictionary with fields for statistics for that gene
    """
    gene_ids_dict = dict()
    gff_data = gpf.ll(annotations_gff_path)
    for line in gff_data:
        if line.startswith("#") == False:
            split_line = line.split()
            if split_line[2] == query_gff_feature:

                gene_id = split_line[8]
                gene_id = gene_id.split("ID=")[1]
                gene_id = gene_id.split(";")[0]
                gene_ids_dict[gene_id] = {"scaff": None, "start": None, "end": None, "strand": None, "orthomcl_cluster": None, "paralog_count": 0, "ortholog_count": 0, "protein_conservation_ratio": 0}

                gene_ids_dict[gene_id]["scaff"] = split_line[0]
                gene_ids_dict[gene_id]["start"] = split_line[3]
                gene_ids_dict[gene_id]["end"] = split_line[4]
                gene_ids_dict[gene_id]["strand"] = split_line[6]
    return gene_ids_dict


def get_orthomcl_species_list(orthomcl_data):
    """
    Gets the list of species identifiers from the OrthoMCL file
    """
    species_ids = list()
    for line in orthomcl_data:
        line = line.split(":")[1]
        line = line.strip()
        split_line = line.split()
        split_line = [n for n in split_line if "(" in n and ")" in n]
        split_line = [n.split("(")[1] for n in split_line]
        split_line = [n.split(")")[0] for n in split_line]
        for item in split_line:
            if item not in species_ids:
                species_ids.append(item)
    return species_ids


def get_species_counts_in_cluster_line(cluster_line, species_ids):
    """
    Input: 1) a line from OrthoMCL output file, 2) list of unique species IDs in the OrthoMCL file
    Output: dictionary where the keys are species IDs and the values are the counts of proteins of the selected species in the selected cluster
    """
    species_counts_dict = dict()
    for item in species_ids:
        species_counts_dict[item] = 0
    cluster_line = cluster_line.split("\t")[1]
    split_cluster_line = cluster_line.split()
    split_cluster_line = [n.split("(")[1] for n in split_cluster_line]
    split_cluster_line = [n.split(")")[0] for n in split_cluster_line]
    for species_id in species_counts_dict:
        species_counts_dict[species_id] = split_cluster_line.count(str(species_id))
    return species_counts_dict


def main(query_species_id, annotations_gff_path, orthomcl_file_path, out_path, query_gff_feature):
    gene_ids_dict = get_gene_ids_dict(annotations_gff_path, query_gff_feature)
    orthomcl_data = gpf.l(orthomcl_file_path)
    species_ids = get_orthomcl_species_list(orthomcl_data)
    nr_of_species = len(species_ids)
    gene_not_found_errors_count = 0

    for line in orthomcl_data:
        cluster_name = line.split("(")[0]
        species_counts_dict = get_species_counts_in_cluster_line(line, species_ids)
        if species_counts_dict[query_species_id] > 0:
            cluster_paralog_count = species_counts_dict[query_species_id] - 1
            cluster_ortholog_count = sum(species_counts_dict.values()) - species_counts_dict[query_species_id]
            conservation_list = [n for n in species_counts_dict.values() if n > 0]
            protein_conservation_ratio = (len(conservation_list) - 1) / (nr_of_species - 1)

            cluster_line = line.split("\t")[1]
            split_cluster_line = cluster_line.split()
            for item in split_cluster_line:
                item_species_id = item.split("(")[1]
                item_species_id = item_species_id.split(")")[0]
                if item_species_id == query_species_id:
                    item_gene_id = item.split("(")[0]
                    if item_gene_id not in gene_ids_dict:
                        sys.stderr.write("Warning: " + item_gene_id + " occurs in the OrthoMCL file but was not found annotated as a gene in the input GFF file\n")
                        gene_not_found_errors_count += 1
                        if gene_not_found_errors_count == 5:
                            sys.stderr.write("Mismatch between gene names in the OrthoMCL file ({}) and input GFF file ({}) appears at least 5 times\n".format(orthomcl_file_path, annotations_gff_path))
                            sys.exit(1)
                    else:
                        gene_ids_dict[item_gene_id]["orthomcl_cluster"] = cluster_name
                        gene_ids_dict[item_gene_id]["paralog_count"] = cluster_paralog_count
                        gene_ids_dict[item_gene_id]["ortholog_count"] = cluster_ortholog_count
                        gene_ids_dict[item_gene_id]["protein_conservation_ratio"] = protein_conservation_ratio
    df = pd.DataFrame(gene_ids_dict)
    df = df.transpose()
    df.index.name = "gene"
    df.reset_index(inplace=True)
    df.to_csv(out_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query_species_id", type=str, help="Query species ID in OrthoMCL file. This is the species for which the output CSV is generated. \
        The species IDs are in brackets after protein IDs in the OrthoMCL file. \
        E.g. if a line in the OrthoMCL contains 'EMAXIMA_000505700.1(emaxima) EMAXIMA_000512600.1(emaxima)', the species ID for these proteins is emaxima")
    parser.add_argument("annotations_gff_path", type=str, help="Path to GFF3 file with gene annotations for the query species")
    parser.add_argument("orthomcl_file_path", type=str, help="Path to OrthoMCL output file that contains the OrthoMCL clusters for the query species and the related species")
    parser.add_argument("out_path", type=str, help="Path for output CSV file")
    parser.add_argument("--query_gff_feature", type=str, help="The tag for gene features in the GFF3 file (default: 'mRNA')", default="mRNA")
    args = parser.parse_args()
    main(args.query_species_id, args.annotations_gff_path, args.orthomcl_file_path, args.out_path, args.query_gff_feature)


