#!/usr/bin/env python3
"""
Script for counting cluster junctions between windows in the GDA output 'clusters.bed' file.
Fisher test is used to determine if some types of cluster junctions occur more rarely or often that expected by chance
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

from collections import OrderedDict
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest
import general_purpose_functions as gpf
import argparse


def get_observed_counts(clusters_dict, cluster_pairs):
    """
    Returns the counts of cluster junctions that can be seen between elements of clusters_list
    """
    observed_counts_dict = OrderedDict()
    for cluster_pair in cluster_pairs:
        observed_counts_dict[(cluster_pair[0], cluster_pair[1])] = 0
    for scaff in clusters_dict:
        scaff_clusters_list = clusters_dict[scaff]
        for pos in range(0, len(scaff_clusters_list)-1):
            pair = (scaff_clusters_list[pos], scaff_clusters_list[pos+1])
            observed_counts_dict[pair] += 1
    return observed_counts_dict


def get_probablities_dict(clusters_list, cluster_pairs):
    """
    Finds the probabilities of observing each type of cluster junctions if the junctions were distributed by chance,
        given the number of windows in each cluster
    """
    clusters_list_len = len(clusters_list)
    probabilities_dict = OrderedDict()
    probabilities_sum = 0
    list_for_prob2_len = clusters_list_len -1
    for cluster_pair in cluster_pairs:
        item1 = cluster_pair[0]
        item2 = cluster_pair[1]
        prob1 = clusters_list.count(item1)/clusters_list_len
        list_for_prob2 = clusters_list.copy()
        list_for_prob2.remove(item1)
        prob2 = list_for_prob2.count(item2)/list_for_prob2_len
        pair_probability = prob1 * prob2
        probabilities_sum += pair_probability

        if cluster_pair in probabilities_dict:
            probabilities_dict[cluster_pair] += pair_probability
        else:
            probabilities_dict[cluster_pair] = pair_probability
    assert round(probabilities_sum, 5) == 1.0 # The probabilities should add up to 1. The rounding is to compensate for the imprecision inherent to floating point numbers
    return probabilities_dict


def get_cluster_pairs(clusters_list):
    """
    Input: list of clusters
    Output: all unique pairs that can be formed between clusters of each type. Self-pairing is allowed
    """
    unique_clusters_list = set(clusters_list)
    cluster_pairs = list()
    for cl1 in unique_clusters_list:
        for cl2 in unique_clusters_list:
            cluster_pairs.append((cl1, cl2))
    return cluster_pairs


def cluster_junctions_fisher_test(observed_counts_dict, expected_counts_dict):
    """
    Does a Fisher test of observed vs expected junction count values. Returns a dataframe with cluster names (index), observed and expected junction counts,
        Fisher test odds ratio and p values

    Fisher test matrix:

                       observed    expected
     selected_junction a           b
     other_junctions   c           d
    """
    df_observed = pd.DataFrame.from_dict(observed_counts_dict, orient='index')
    df_expected = pd.DataFrame.from_dict(expected_counts_dict, orient='index')
    df_observed.columns = ["observed"]
    df_expected.columns = ["expected"]

    df = pd.merge(df_observed, df_expected, left_index=True, right_index=True)
    df["expected_int"] = round(df["expected"], 0).astype(int)
    observed_sum = df["observed"].sum()
    expected_int_sum = df["expected_int"].sum()

    oddsratios_list = list()
    pvalues_list = list()
    for row_tuple in df.iterrows():
        row = row_tuple[1]
        a = int(row["observed"])
        b = int(row["expected_int"])
        c = observed_sum - a
        d = expected_int_sum - b
        fisher_matrix = [[a, b], [c, d]]
        oddsratio, pvalue = stats.fisher_exact(fisher_matrix)
        oddsratios_list.append(oddsratio)
        pvalues_list.append(pvalue)

    corrected_pvalues = statsmodels.stats.multitest.multipletests(pvalues_list, alpha=0.05, method="fdr_bh")

    df["oddsratio"] = oddsratios_list
    df["pvalue"] = pvalues_list
    df["corrected_pvalue"] = corrected_pvalues[1]
    df["significance_bool"] = corrected_pvalues[0]
    return df


def load_clusters_dict(bed_file_path):
    """
    Loads cluster positions information from BED file
    """
    clusters_dict = OrderedDict()
    bed_data = gpf.l(bed_file_path)
    for line in bed_data:
        split_line = line.split()
        scaff = split_line[0]
        if scaff not in clusters_dict:
            clusters_dict[scaff] = list()
        cluster = int(split_line[3])
        clusters_dict[scaff].append(cluster)
    return clusters_dict


def condense_counts_dict(input_dict, cluster_pairs_condensed):
    """
    Goes through the dictionary with cluster junction counts and merges the counts for cluster pairs that have the same elements.
        E.g. pair (1,4) is the same as (4,1), so the counts of (4,1) added to the counts of (1,4), and (4,1) is removed from the dictionary
    """
    condensed_dict = OrderedDict()
    for item in input_dict:
        if item in cluster_pairs_condensed:
            condensed_dict[item] = input_dict[item]
        else:
            reversed_item = (item[1], item[0])
            condensed_dict[reversed_item] += input_dict[item]
    return condensed_dict


def main(bed_file_path, out_folder):
    gpf.run_system_command("mkdir -p {}".format(out_folder))
    clusters_dict = load_clusters_dict(bed_file_path)

    clusters_list = list(clusters_dict.values())
    clusters_list = [item for sublist in clusters_list for item in sublist] # Flattening a list of lists

    cluster_pairs = get_cluster_pairs(clusters_list)

    probabilities_dict = get_probablities_dict(clusters_list, cluster_pairs)
    observed_counts_dict = get_observed_counts(clusters_dict, cluster_pairs)
    observed_counts_sum = sum(observed_counts_dict.values())
    expected_counts_dict = OrderedDict()
    for item in observed_counts_dict:
        expected_counts_dict[item] = probabilities_dict[item] * observed_counts_sum

    cluster_pairs_condensed = list()
    for item in cluster_pairs:
        reversed_item = (item[1], item[0])
        if reversed_item not in cluster_pairs_condensed:
            cluster_pairs_condensed.append(item)

    observed_counts_dict2 = condense_counts_dict(observed_counts_dict, cluster_pairs_condensed)
    expected_counts_dict2 = condense_counts_dict(expected_counts_dict, cluster_pairs_condensed)

    fisher_df = cluster_junctions_fisher_test(observed_counts_dict2, expected_counts_dict2)
    fisher_df["observed_vs_expected_fold_diff"] = fisher_df["observed"] / fisher_df["expected_int"]
    fisher_df["log2_observed_vs_expected_fold_diff"] = np.log2(fisher_df["observed_vs_expected_fold_diff"])

    fisher_out_path = out_folder + "/cluster_junctions_fisher_test.csv"
    fisher_df.to_csv(fisher_out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_file_path", type=str, help="Path to the input clusters.bed file")
    parser.add_argument("out_folder", type=str, help="Folder path for output CSV file")
    args = parser.parse_args()
    main(args.bed_file_path, args.out_folder)
