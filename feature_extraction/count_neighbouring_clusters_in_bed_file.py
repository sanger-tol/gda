#!/usr/bin/env python3
"""
Script for counting neighboring clusters of each cluster in GDA output BED file
"""

import general_purpose_functions as gpf
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def get_max_cluster_nr(in_data):
    """
    Returns the largest cluster number in the input BED file
    """
    max_cluster = -1
    for line in in_data:
        split_line = line.split()
        cluster = int(split_line[3])
        if cluster > max_cluster:
            max_cluster = cluster
    return max_cluster


def main(in_path, out_path):

    in_data = gpf.l(in_path)
    max_cluster = get_max_cluster_nr(in_data)
    clusters_array = np.zeros((max_cluster + 2, max_cluster +2), dtype=int)
    df = pd.DataFrame(clusters_array)

    df.columns = (range(-1, max_cluster + 1))
    df.index = (range(-1, max_cluster + 1))

    previous_scaff = None
    previous_cluster = None

    for line in in_data:
        
        split_line = line.split()
        cluster = int(split_line[3])
        scaff = split_line[0]
        if scaff == previous_scaff:
            if previous_cluster is not None:
                if cluster != previous_cluster:
                    
                    junction = None
                    if previous_cluster < cluster:
                        junction = (previous_cluster, cluster)
                    else:
                        junction = (cluster, previous_cluster)
                    df.loc[junction[0], junction[1]] += 1
                    df.loc[junction[1], junction[0]] += 1
        previous_scaff = scaff
        previous_cluster = cluster


    df = (df.apply(np.log2))

    cmap=LinearSegmentedColormap.from_list('rg',["r", "y", "g"], N=256) 

    fig = plt.figure()

    fig.suptitle("Junctions between UMAP clusters in BED file", fontsize=15)
    plt.pcolor(df, cmap=cmap)
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel("log2(number of junctions)")
    plt.savefig(out_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to input BED file that has been produced by the clustering script of GDA")
    parser.add_argument("out_path", type=str, help="Path for output plot PNG file")
    args = parser.parse_args()
    main(args.in_path, args.out_path)
