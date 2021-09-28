#!/usr/bin/env python3
"""
Script for processing a GFF file to add missing IDs. The input is a GFF that has been created by combining Augustus, Barrnap and tRNAscan-SE
"""

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf

in_path = sys.argv[1]
species_id = sys.argv[2]

in_data = gpf.ll(in_path)
counters_dict = {"exon": 1, "intron": 1, "rRNA": 1, "start_codon": 1, "stop_codon": 1}

for line in in_data:
    if line != "###":
        if line.startswith("#") == False:
            split_line = line.split()
            feature_type = split_line[2]
            split_line[8] = split_line[8].replace("ID=", "ID=" + species_id + "_")
            split_line[8] = split_line[8].replace("Parent=", "Parent=" + species_id + "_")
            if feature_type == "CDS":
                print("\t".join(split_line))
                exon_split_line = split_line
                exon_name = "exon-" + str(counters_dict["exon"])
                exon_parent = exon_split_line[8].split(";")[1]
                exon_split_line[7] = "."
                exon_split_line[8] = "ID=" + exon_name + ";" + exon_parent
                exon_split_line[2] = "exon"
                print("\t".join(exon_split_line))
                counters_dict["exon"] += 1
            elif feature_type == "tRNA":
                gene_split_line = split_line.copy()
                gene_split_line[2] = "gene"
                trna_gene_name = gene_split_line[8].split("ID=")[1]
                gene_split_line[8] = "ID=gene-" + trna_gene_name
                split_line[8] = "ID=trna-" + trna_gene_name + ";Parent=gene-" + trna_gene_name
                print("\t".join(gene_split_line))
                print("\t".join(split_line))
            elif feature_type == "rRNA":
                gene_split_line = split_line.copy()
                gene_split_line[2] = "gene"
                rrna_gene_name = feature_type.lower() + "-" + str(counters_dict[feature_type])
                gene_split_line[8] = "ID=gene-" + rrna_gene_name + ";" + gene_split_line[8]
                split_line[8] = "ID=rrna-" + rrna_gene_name + ";Parent=gene-" + rrna_gene_name + ";" + split_line[8]
                print("\t".join(gene_split_line))
                print("\t".join(split_line))
                counters_dict[feature_type] += 1
            elif feature_type in ("intron", "start_codon", "stop_codon"):
                split_line[8] = "ID=" + feature_type.lower() + "-" + str(counters_dict[feature_type]) + ";" + split_line[8]
                print("\t".join(split_line))
                counters_dict[feature_type] += 1
            else:
                print("\t".join(split_line))
        else:
            print(line)

