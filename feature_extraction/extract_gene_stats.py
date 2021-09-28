#!/usr/bin/env python3
"""
Script for extracting gene statistics (length, exon count, intron count, strand) from a GFF file
"""

import numpy as np
import sys
import itertools
import general_purpose_functions as gpf
from collections import OrderedDict
import argparse


class Gene:
    """
    Class that contains the information on the exons of one gene
    """
    def __init__(self, name, scaff, start, end, strand):
        if start >= end:
            raise Exception("Error in gene {} coordinates: start coordinate is larger than end coordinate".format(name))
        if strand == "+":
            strand = 1
        elif strand == "-":
            strand = -1
        self.name = name
        self.scaff = scaff
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = list()
    
    def add_exon(self, exon_start, exon_end):
        """
        Adding one exon to the gene
        """
        if exon_start >= self.start and exon_end <= self.end and exon_start <= exon_end:
            self.exons.append((exon_start, exon_end))
        else:
            sys.stderr.write("Error adding exon to gene {}: exon coordinates do not match the mRNA coordinates\n".format(self.name))
            sys.stderr.write("mRNA start: {}, mRNA end: {}. Exon start: {}, exon end: {}\n".format(str(self.start), str(self.end), str(exon_start), str(exon_end)))


    def get_exon_lengths(self):
        """
        Returns exon lengths as a list
        """
        exon_lengths = [n[1] - n[0] + 1 for n in self.exons]
        return exon_lengths

    def get_intron_lengths(self):
        """
        Returns intron lenghts as a list
        """
        gene_length = self.end + 1 - self.start
        gene_array = np.zeros(gene_length, dtype=bool)
        for exon in self.exons:
            exon_start = exon[0] - self.start + 1
            exon_end = exon[1] - self.start + 1
            gene_array[exon_start - 1:exon_end] = 1
        intron_lengths = [sum(1 for _ in group) for key, group in itertools.groupby(gene_array) if key == False] # https://stackoverflow.com/questions/24342047/count-consecutive-occurences-of-values-varying-in-length-in-a-numpy-array
        return intron_lengths


def print_genes_dict_contents(genes_dict):
    """
    Goes through the genes_dict dictionary that contains the gene statistics and prints them out as a comma separated table
    """
    print("gene,scaff,start,end,strand,gene_length,exon_count,gene_average_exon_length,gene_average_intron_length")
    for mrna_id in genes_dict:
        mrna_gene = genes_dict[mrna_id]
        gene_exon_count = len(mrna_gene.exons)
        if gene_exon_count > 0:
            gene_average_exon_length = np.mean(mrna_gene.get_exon_lengths())
            intron_lengths = mrna_gene.get_intron_lengths()
            gene_average_intron_length = np.nan
            if len(intron_lengths) != 0:
                gene_average_intron_length = np.mean(mrna_gene.get_intron_lengths())
            gene_length = mrna_gene.end - mrna_gene.start + 1
            out_list = [mrna_gene.name, mrna_gene.scaff, str(mrna_gene.start), str(mrna_gene.end), str(mrna_gene.strand), str(gene_length), str(gene_exon_count), str(gene_average_exon_length), str(gene_average_intron_length)]
            print(",".join(out_list))
        else:
            sys.stderr.write("No exons with valid coordinates found for mRNA {}\n".format(str(mrna_id)))


def main(gff_path):    
    genes_dict = OrderedDict()
    gff_data = gpf.l(gff_path)
    for line in gff_data:
        if line.startswith("#") == False:
            split_line = line.split()
            field8 = split_line[8]
            if split_line[2] == "mRNA":
                mrna_id = gpf.spl(field8, "ID=", ";")
                mrna_gene = Gene(mrna_id, split_line[0], int(split_line[3]), int(split_line[4]), split_line[6])
                genes_dict[mrna_id] = mrna_gene
            elif split_line[2] == "exon":
                exon_gene_id = gpf.spl(field8, "Parent=", ";")
                exon_start = int(split_line[3])
                exon_end = int(split_line[4])
                if exon_gene_id in genes_dict:
                    genes_dict[exon_gene_id].add_exon(exon_start, exon_end)
    print_genes_dict_contents(genes_dict)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff_path", type=str, help="Path to input GFF file")
    args = parser.parse_args()
    main(args.gff_path)
    