#!/usr/bin/env python3
"""
This script will iterate over N-neighbours and cluster_size parameters for the dimension reduction and clustering steps of GDA.
It will report the resulting clusters allowing the user to select parameters for a full run of the analysis step.
"""
# MIT License
# 
# Copyright (c) 2020-2021 Genome Research Ltd.
# 
# Authors: Adam Reid (ar11@sanger.ac.uk), Eerik Aunin (ea10@sanger.ac.uk)
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

import argparse
import sys
import os
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import umap
import hdbscan
import matplotlib.pyplot as plt
import spectra
import numpy as np
from matplotlib.colors import ListedColormap

from gda_clustering import get_palette, run_umap, run_hdbscan, check_umap_version
from sklearn import metrics

# Ignore NUMBA warnings related to running UMAP
import warnings
from numba.core.errors import NumbaPerformanceWarning
from numba.core.errors import NumbaDeprecationWarning
from numba.core.errors import NumbaWarning
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings("ignore", category=NumbaWarning)
# Ignore matplotlib's "More than 20 figures have been opened" warning
plt.rcParams.update({'figure.max_open_warning': 0})



###############
# Functions
###############


def get_cluster_props(cluster_list):
    '''Makes a dict of cluster number to genome proportion'''
    cluster_genome_prop = dict()
    cluster_genome_prop[-1] = 0
    for x in sorted(set(cluster_list)):
        cluster_size = list(cluster_list).count(x)
        cluster_genome_prop[x] = (cluster_size * 100) / cluster_list.size
    return(cluster_genome_prop)


def print_clustering_result(n, c, cluster_labels, metrics_dict_entry):
    '''Prints clustering result to STDOUT'''
    print('N neighbours: {}\nMin cluster size: {}'.format(n, c))
    for x in cluster_labels:
        print('Cluster {} is {:.2f}% of the genome'.format(x, cluster_labels[x]))

    silhouette_score = metrics_dict_entry["silhouette_score"]
    if silhouette_score is not None:
        print("Silhouette score: {:.2f}".format(silhouette_score))
        print("Davies-Bouldin index: {:.2f}".format(metrics_dict_entry["davies_bouldin_index"]))
        print("Calinski-Harabasz score: {:.2f}".format(metrics_dict_entry["calinski_harabasz_score"]))
    print()


def set_up_outdir(outdir):
    '''Sets up output directory'''
    if not os.path.exists(outdir):
        sys.stderr.write("Creating directory: {}\n".format(outdir))
        os.mkdir(outdir)
    else:
        sys.stderr.write("{} already exists\n".format(outdir))

    if not os.path.exists(outdir + '/parameter_selection'):
        sys.stderr.write("Creating directory: {}\n".format(outdir + '/parameter_selection'))
        os.mkdir(outdir + '/parameter_selection')
    else:
        sys.stderr.write("{} already exists\n".format(outdir + '/parameter_selection'))


def get_html_string(neighbours_list, cluster_size_list, plot_files, all_cluster_props, silhouette_scores):
    '''Produces a string of HTML code that will later be written as a file'''
    html = '<html><body>\n'
    html = html + '<font face="Arial">\n'

    # Make a 2x2 table with table of figures in last element
    # First row
    html = html + '<table border=0><tr><td></td><td style="text-align:center"><h1>N neighbours</h1></td></tr>\n'

    # Second row
    html = html + '<tr><td><h1>Min cluster size</h1></td><td>'

    # Table of figures
    html = html + '<table border=0><tr><th></th>'
    for n in neighbours_list:
        html = html + '<th><h2>{}</h2></th>'.format(n)
    html = html + '</tr>\n'

    workdir = os.getcwd()

    for c in cluster_size_list:
        html = html + '<tr><td><h2>{}</h2></td>'.format(c)
        for n in neighbours_list:
            img_file_basename = os.path.basename(workdir + "/" + plot_files[c][n])
            silhouette_score = silhouette_scores[n][c]
            if silhouette_score is not None:
                silhouette_score = round(silhouette_score, 2)
            silhouette_score = str(silhouette_score)
            html = html + '<td><table><tr><td><img src="{}"></td></tr><tr><td style="text-align:center"><h2>{:.2f}% unclassified<br>Silhouette score: {}</h2></td></tr></table></td>'.format(img_file_basename, all_cluster_props[n][c][-1], silhouette_score)
        html = html + '</tr>\n'

    html = html + '</table>'

    html = html + " </font>\n"
    # End enclosing table
    html = html + '</tr></table></body></html>'
    return html


def get_neighbours_list(n):
    '''Extracts and check n_neighbours values from a command line argument'''
    neighbours_list = n.split(',')
    # Check whether values are all numeric
    test = [i.isnumeric() for i in neighbours_list]
    if test.count(True) != len(neighbours_list):
        sys.stderr.write('Not all N neighbours values are numeric: {}\n'.format(neighbours_list))
        sys.exit(1)
    else:
        neighbours_list = [int(x) for x in neighbours_list]
    return neighbours_list


def get_cluster_size_list(c):
    '''Extracts and checks cluster size values from a command line argument'''
    cluster_size_list = c.split(',')

    # Check whether values are all numeric
    test = [i.isnumeric() for i in cluster_size_list]
    if test.count(True) != len(cluster_size_list):
        sys.stderr.write('Not all cluster size values are numeric: {}\n'.format(cluster_size_list))
        sys.exit(1)
    else:
        cluster_size_list = [int(x) for x in cluster_size_list]
    return cluster_size_list


def filter_df_to_keep_selected_scaff(selected_scaff, data):
    '''If selected_scaff is not None, filters the data from the TSV file to keep only the entries for one selected scaffold'''
    if selected_scaff is not None:
        if selected_scaff in data['chromosome'].tolist():
            data = data[data['chromosome'] == selected_scaff]
        else:
            sys.stderr.write("Scaffold with the name '{}' was not found\n".format(selected_scaff))
            sys.exit(1)
    return data


def export_metrics_table(metrics_dict, outdir_full):
    '''Exports the clustering metrics as a CSV table'''
    metrics_df = pd.DataFrame(metrics_dict)
    metrics_df = metrics_df.transpose()
    metrics_df = metrics_df.sort_values(["silhouette_score", "unclassified_percentage", "davies_bouldin_index", "calinski_harabasz_score"], ascending=[False, True, True, False])
    metrics_outfile = outdir_full + "/clustering_metrics.csv"
    metrics_df.to_csv(metrics_outfile, index=True)
    


def main(args):
    #################
    # Procedural code
    #################
    
    check_umap_version()

    neighbours_list = get_neighbours_list(args.n_neighbors)
    cluster_size_list = get_cluster_size_list(args.cluster_size_cutoff)
    outdir = args.directory
    tracks_file = args.tracks

    set_up_outdir(outdir)
    outdir_full = outdir + '/parameter_selection'

    # Read in GDA data tracks as pandas data frame
    data = pd.read_csv(tracks_file, sep='\t', index_col=0)
    data = filter_df_to_keep_selected_scaff(args.selected_scaff_only, data)

    # Filter windows remove any N-containing windows
    if 'N_percentage' in data.columns:
        data.drop(index=data[data.N_percentage>0].index, axis=0, inplace=True)

    # Make a dataframe with only window_name index and feature values suitable for clustering
    data_to_cluster = data.drop(['species', 'start', 'end', 'chromosome'], axis=1)

    plot_files = dict()


    # Save cluster props
    all_cluster_props = dict()

    silhouette_scores = dict()

    metrics_dict = dict()

    # Loop over n_neighbours values and run UMAP
    for n in neighbours_list:

        # Get UMAP embedding for the data
        embedding = run_umap(data_to_cluster, n)

        # Loop over cluster window sizes and run HDBSCAN
        for c in cluster_size_list:
            cluster_labels = run_hdbscan(embedding, args.leaf_size, args.min_samples, min_cluster_size=c)

            # Calculate Silhouette score, Davies-Bouldin index and Calinski-Harabasz score to assess how well clustering with the selected settings has worked
            silhouette_score = None
            davies_bouldin_index = None
            calinski_harabasz_score = None

            if len(set(cluster_labels)) > 1:
                silhouette_score = metrics.silhouette_score(embedding, cluster_labels)
                davies_bouldin_index = metrics.davies_bouldin_score(embedding, cluster_labels)
                calinski_harabasz_score = metrics.calinski_harabasz_score(embedding, cluster_labels)

            cluster_props = get_cluster_props(cluster_labels)

            if n not in all_cluster_props:
                all_cluster_props[n] = dict()
            all_cluster_props[n][c] = cluster_props

            if n not in silhouette_scores:
                silhouette_scores[n] = dict()
            silhouette_scores[n][c] = silhouette_score

            metrics_dict_entry = {"n_neighbors": n, "min_cluster_size": c, "silhouette_score": silhouette_score, "unclassified_percentage": cluster_props[-1], "davies_bouldin_index": davies_bouldin_index, "calinski_harabasz_score": calinski_harabasz_score}
            metrics_dict_entry_id = str(n) + "_" + str(c)
            metrics_dict[metrics_dict_entry_id] = metrics_dict_entry

            print_clustering_result(n, c, cluster_props, metrics_dict_entry)

            palette = get_palette(max(cluster_labels) + 1)
            if -1 in cluster_labels:
                palette.insert(0, "#B3B4B5")
            cmap1 = ListedColormap(palette)

            ax = plt.subplots()[1]
            scatter = ax.scatter(embedding[:,0], embedding[:,1], c=cluster_labels, cmap=cmap1, s=0.3)

            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            legend1 = ax.legend(*scatter.legend_elements(),
                        loc="center left", title="Clusters", bbox_to_anchor=(1, 0.5))
            ax.add_artist(legend1)

            outfile = '{}/umap_clustering_n{}_c{}.png'.format(outdir_full, n, c)
            plt.savefig(outfile)
            plot_file = outfile

            if c not in plot_files:
                plot_files[c] = dict()
            plot_files[c][n] = plot_file


    html = get_html_string(neighbours_list, cluster_size_list, plot_files, all_cluster_props, silhouette_scores)

    with open('{}/parameters.html'.format(outdir_full), 'w') as o:
        o.write(html)

    sys.stderr.write('Please load this file into your browser to view results:\n\n{}/parameters.html\n'.format(outdir_full))

    export_metrics_table(metrics_dict, outdir_full)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome Decomposition Analysis of windowed tracks - parameter selection')
    parser.add_argument("-n", "--n_neighbors", help="N neighbours argument for UMAP [5,10,15,20,50,100]", default='5,10,15,20,50,100', type=str)
    parser.add_argument("-c", "--cluster_size_cutoff", help="HDBSCAN min cluster size [50,100,200,500]", default='50,100,200,500', type=str)
    parser.add_argument("-d", "--directory", help="Output dir [gda_out]", default="gda_out", type=str)
    parser.add_argument("--leaf_size", help="leaf_size setting for HDBSCAN (default: 40)", default=40, type=int)
    parser.add_argument("--min_samples", help="min_samples setting for HDBSCAN (default: None)", default=None, type=int)
    parser.add_argument("--selected_scaff_only", help="Optional: name of a scaffold from the input TSV file, to run UMAP+HDBSCAN with that scaffold only. Default: None (all scaffold are used by default)", default=None, type=str)
    parser.add_argument("tracks", help="File of windowed GDA tracks", type=str)
    args = parser.parse_args()
    main(args)
