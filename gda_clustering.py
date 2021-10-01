#!/usr/bin/env python3
"""
Code to generate GDA analysis
A script for turning GDA tracks into analysed clusters
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

# Several species-specific analyses are done based on BedFile objects
from BedFile import BedFile

import sys
import os
import argparse
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import MinMaxScaler
import hdbscan
from scipy.stats import ks_2samp
import json
import spectra

# Ignore NUMBA warnings related to running UMAP
import warnings
import numba
from numba import NumbaPerformanceWarning
from numba import NumbaDeprecationWarning
from numba import NumbaWarning
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings("ignore", category=NumbaWarning)

###############
# Functions
###############

def get_palette(n):
    '''Returns a palette equally spaced hues with n colours'''
    hues = np.linspace(15, 375, n+1).tolist()
    hues = hues[0:len(hues)-1]
    cols = [spectra.lch(65, 100, n) for n in hues]
    hexcols = [n.hexcode for n in cols]
    return hexcols


def get_cluster_cols_dict(max_n):
    '''Returns a dictionary where the keys are cluster IDs and the values are colours assigned to the clusters'''
    palette = get_palette(max_n + 1)
    palette.insert(0, "#B3B4B5")
    colours_dict = dict()
    counter = 0
    for i in range(-1, max_n + 1):
        colours_dict[i] = palette[counter]
        counter += 1
    return colours_dict


def run_umap(data, n_neighbors=13):
    '''Scale feature values each from 0 to 1 and run UMAP, return 2D embedding'''
    # Scale data 0->1
    scalar = MinMaxScaler()
    scalar.fit(data)
    scaled_data = scalar.transform(data)
    reducer = umap.UMAP(random_state=123, n_neighbors=n_neighbors, min_dist=0.1, n_components=2, metric='euclidean')
    embedding = reducer.fit_transform(scaled_data)
    return(embedding)


def run_hdbscan(embedding, leaf_size, min_samples, min_cluster_size=200):
    '''Cluster UMAP embedding with HDBSCAN, return cluster labels'''
    clusterer = hdbscan.HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True,
        gen_min_span_tree=False, leaf_size=leaf_size,
        metric='euclidean', min_cluster_size=min_cluster_size, min_samples=min_samples, p=None)
    clusterer.fit(embedding)
    return(clusterer.labels_)


def write_cluster_bed(data, cluster_colours, outdir):
    '''Write BED file describing clusters - return dict of species -> filename'''
    # Declare list of filnames
    clusters_bed_fn = ''
    # Loop over species to make BED file for each
    clusters_bed_fn = ('clusters.bed')
    cb = open(outdir + '/' + clusters_bed_fn, 'w')
    for row_tuple in data.iterrows():
        row = row_tuple[1]
        bed_string = '\t'.join([row.chromosome, str(row.start), str(row.end + 1), str(row.cluster), '0', '+', str(row.start), str(row.end + 1), cluster_colours[row.cluster]])
        cb.write(bed_string + "\n")
    cb.close()
    return(clusters_bed_fn)


def write_cluster_gff(data, outdir):
    '''Write GFF of clusters'''
    clusters_gff_fn = 'clusters.gff'
    cg = open(outdir + '/' + clusters_gff_fn, 'w')
    for row_tuple in data.iterrows():
        row = row_tuple[1]
        start = row.start + 1
        end = row.end + 1
        gff_string = '\t'.join([row.chromosome, 'GDA', 'region', str(start), str(end), '.', '+', '.', 'ID='+str(row.cluster)+';colour=1'])
        cg.write(gff_string + "\n")
    cg.close()


def get_significant_features(data, pvalue_cutoff):
    '''Get significant features in each cluster and determine proportions of each cluster in each species'''
    cluster_list = sorted(set(data['cluster'].tolist()))

    # Make a dictionary of the data for ease of access
    data_dict = data.to_dict()
    cluster_dict = data_dict['cluster']

    # Make a dictionary with cluster as key and a list of windows as the value
    cluster_to_windows = dict()
    for x in cluster_dict:
        if cluster_dict[x] not in cluster_to_windows:
            cluster_to_windows[cluster_dict[x]] = list()
        cluster_to_windows[cluster_dict[x]].append(x)

    # Record a set of significant features for returning
    sig_features = set()
    cluster_means = dict() # This is used for heatmap
    cluster_genome_prop = dict() # Proportion of genome made up by cluster
    # cluster proportions by species
    cluster_species_genome_prop = dict()

    species_list = sorted(set(data['species'].tolist()))

    # Store signficant feature results for each cluster
    cluster_results = dict()
    cluster_results_df = pd.DataFrame(columns=['cluster', 'feature', 'cluster_data.size', 'other_data.size', 'stat_less', 'pvalue_less', 'stat_great', 'pvalue_great', 'cluster_median', 'other_median', 'cluster_mean', 'other_mean'])

    for c in cluster_list:
        cluster_data = data[data['cluster'] == c]
        # Store proportion of genome made up of each cluster
        cluster_genome_prop[c] = (cluster_data.size * 100) / (data.size)

        # Determine proportion of cluster in each species' genomes
        for s in species_list:
            species_cluster_data = data[data.cluster.eq(c) & data.species.eq(s)]
            if c not in cluster_species_genome_prop:
                cluster_species_genome_prop[c] = dict()
            cluster_species_genome_prop[c][s] = (species_cluster_data.size * 100) / (data[data.species.eq(s)].size)

        for f in data.columns:
            if f in ['species', 'cluster', 'start', 'end', 'chromosome']:
                continue
            # select values for this feature and this cluster
            cluster_data = data[f][cluster_to_windows[c]]
            other_index = data.index.difference(cluster_to_windows[c])
            other_data = data[f][other_index]

            # Two-sided kstest seems to give a lot of pvalue=1 where it shouldn't
            (stat_less, pvalue_less) = ks_2samp(other_data, cluster_data, alternative="less", mode="exact")
            (stat_great, pvalue_great) = ks_2samp(other_data, cluster_data, alternative="greater", mode="exact")
            cluster_median = np.median(cluster_data)
            other_median = np.median(other_data)
            cluster_mean = np.mean(cluster_data)
            other_mean = np.mean(other_data)

            #Store cluster means for heatmap
            if f not in cluster_means:
                cluster_means[f] = dict()
            cluster_means[f][c] = cluster_mean

            if pvalue_less <= pvalue_cutoff or pvalue_great <= pvalue_cutoff:
                sig_features.add(f)
                if c not in cluster_results:
                     cluster_results[c] = dict()
                if f not in cluster_results[c]:
                     cluster_results[c][f] = list()
                cluster_results[c][f] = [cluster_data.size, other_data.size, stat_less, pvalue_less, stat_great, pvalue_great, cluster_median, other_median, cluster_mean, other_mean]
                cluster_results_df = cluster_results_df.append({'cluster' : c, 'feature' : f, 'cluster_data.size' : cluster_data.size, 'other_data.size' : other_data.size, 'stat_less' : '{:.2f}'.format(stat_less), 'pvalue_less' : '{:.2e}'.format(pvalue_less), 'stat_great' : '{:.2f}'.format(stat_great), 'pvalue_great' : '{:.2e}'.format(pvalue_great), 'cluster_median' : '{:.5f}'.format(cluster_median), 'other_median' : '{:.5f}'.format(other_median), 'cluster_mean' : '{:.5f}'.format(cluster_mean), 'other_mean' : '{:.5f}'.format(other_mean)}, ignore_index=True)
    return(sig_features, cluster_results_df, cluster_species_genome_prop, cluster_means)


def feature_histograms(data, sig_features):
    '''Generate data for feature histograms'''
    feat_hist_dict = dict()
    for f in sig_features:
        for c in sorted(set(data['cluster'].tolist())):
            # Select rows matching feature and cluster from data
            data_subset = data.loc[data.cluster==c]
            if f not in feat_hist_dict:
                feat_hist_dict[f] = dict()
            if c not in feat_hist_dict[f]:
                feat_hist_dict[f][c] = list()
            feat_hist_dict[f][c] = data_subset[f].tolist()
    return(feat_hist_dict)


# Functions for writing data out for display with Shiny
def write_feature_table(table_df, outfile, outdir):
    '''Write table of significant features for reading by Shiny'''
    table_df.to_csv(outdir + '/' + outfile, index=False)


def write_embedding(emb, cluster_labels, species, outfile, outdir, cluster_colours):
    '''Write out UMAP embedding for display with Shiny'''
    o = open(outdir + '/' + outfile, 'w')
    o.write(','.join(['UMAP1', 'UMAP2', 'species', 'cluster', 'color']))
    o.write('\n')
    for i in range(0,len(cluster_labels)):
        o.write(','.join([str(emb[i,0]), str(emb[i,1]), str(species[i]), str(cluster_labels[i]), str(cluster_colours[cluster_labels[i]])]))
        o.write('\n')
    o.close()


def write_cluster_heatmap(cluster_means, sig_features, outfile, outdir):
    '''Produce an output file with the data for the heatmap of genomic features per cluster'''
    norm_vals = dict()
    for f in cluster_means:
        if f in sig_features:
            min_val = min(cluster_means[f].values())
            max_val = max(cluster_means[f].values())
            for c in cluster_means[f]:
                if f not in norm_vals:
                    norm_vals[f] = dict()
                norm_vals[f][c] = (cluster_means[f][c] - min_val) / (max_val - min_val)
    df = pd.DataFrame.from_dict(norm_vals)
    df.to_csv(outdir + '/' + outfile)


def write_genome_prop(genome_props, outfile, outdir):
    '''Write out CSV of proportions of genome composed of each cluster'''
    df = pd.DataFrame.from_dict(genome_props)
    df.to_csv(outdir + '/' + outfile)


def write_chr_comp_heatmap(data, outfile, outdir):
    '''Write out chromosome composition heatmap data as CSV'''
    comp_res_df = pd.DataFrame(data)
    comp_res_df = comp_res_df.fillna(0)
    comp_res_df = comp_res_df.transpose()
    comp_res_df.to_csv(outdir + '/' + outfile)


def write_json(data, outfile, outdir):
    '''Write out JSON file of data'''
    with open(outdir + '/' + outfile, 'w') as file:
        file.write(json.dumps(data))


def write_circos_json(json_dict, outfile, outdir):
    '''Write out a JSON file with chromosomal locations of clusters. This can be used to make a Circos plot or a raster plot'''
    with open(outdir + '/' + outfile, 'w') as file:
        file.write(json.dumps(json_dict))


def validate_clustering_params(param_values_dict):
    '''Function for checking if the clustering parameters selected by the user fit into expected ranges. Prints a warning or exits with an error message if they do not fit'''
    param_limits_dict = {"umap_n_neighbors": {"min": 1, "recommended_min": 2, "recommended_max": 200, "max": None},
    "hdbscan_min_cluster_size": {"min": 1, "recommended_min": None, "recommended_max": None, "max": None},
    "pvalue_cutoff": {"min": 0, "recommended_min": None, "recommended_max": 0.05, "max": 1},
    "cluster_position_histogram_window_number": {"min": 1, "recommended_min": 2, "recommended_max": 100, "max": None},
    "leaf_size": {"min": 1, "recommended_min": 1, "recommended_max": None, "max": None},
    "min_samples": {"min": 1, "recommended_min": 1, "recommended_max": None, "max": None}}
    for param_name, param_value in param_values_dict.items():
        if param_value is not None:
            param_limits = param_limits_dict[param_name]
            if param_limits["min"] is not None:
                if param_value < param_limits["min"]:
                    sys.stderr.write("The value for the parameter {} ({}) should be higher than {}\n".format(param_name, param_value, param_limits["min"]))
                    sys.exit(1)
            if param_limits["recommended_min"] is not None:
                if param_value < param_limits["recommended_min"]:
                    sys.stderr.write("Warning: the value for the parameter {} ({}) is lower than the recommended minimum ({})\n".format(param_name, param_value, param_limits["recommended_min"]))
            if param_limits["max"] is not None:
                if param_value > param_limits["max"]:
                    sys.stderr.write("The value for the parameter {} ({}) should be lower than {}\n".format(param_name, param_value, param_limits["max"]))
                    sys.exit(1)
            if param_limits["recommended_max"] is not None:
                if param_value > param_limits["recommended_max"]:
                    sys.stderr.write("Warning: the value for the parameter {} ({}) is higher than the recommended maximum ({})\n".format(param_name, param_value, param_limits["recommended_max"]))


def read_clustering_metrics_from_csv(clustering_metrics_file):
    '''Reads the recommended clustering settings from the clustering_metrics.csv file that has been generated by gda_parameters.py'''
    selected_n_neighbors = None
    selected_min_cluster_size = None

    if os.path.exists(clustering_metrics_file) == False:
        sys.stderr.write("File not found ({})\n".format(clustering_metrics_file))
        sys.exit(1)

    metrics_df = pd.read_csv(clustering_metrics_file, index_col=0)

    if "n_neighbors" and "min_cluster_size" in metrics_df:
        if metrics_df.shape[0] >= 1:
            selected_n_neighbors = int(metrics_df["n_neighbors"][0])
            selected_min_cluster_size = int(metrics_df["min_cluster_size"][0])
        else:
            sys.stderr.write("The clustering metrics file ({}) appears to contain no data rows\n".format(clustering_metrics_file))
            sys.exit(1)
    else:
        sys.stderr.write("The clustering metrics file ({}) does not have the expected columns\n".format(clustering_metrics_file))
        sys.exit(1)
    return selected_n_neighbors, selected_min_cluster_size


def main(umap_n_neighbors, hdbscan_min_cluster_size, pvalue_cutoff, cluster_position_histogram_window_number, leaf_size, min_samples, outdir, metrics_file_path, tracks_file):
    #################
    # Procedural code
    #################

    if metrics_file_path != "":
        umap_n_neighbors, hdbscan_min_cluster_size = read_clustering_metrics_from_csv(metrics_file_path)

    print(">>>")
    print("umap_n_neighbors", umap_n_neighbors, "hdbscan_min_cluster_size", hdbscan_min_cluster_size)

    param_values_dict = {"umap_n_neighbors": umap_n_neighbors, "hdbscan_min_cluster_size": hdbscan_min_cluster_size, "pvalue_cutoff": pvalue_cutoff, "cluster_position_histogram_window_number": cluster_position_histogram_window_number, "leaf_size": leaf_size, "min_samples": min_samples}
    validate_clustering_params(param_values_dict)

    sys.stderr.write('\t'.join([str(umap_n_neighbors), str(hdbscan_min_cluster_size), str(pvalue_cutoff), tracks_file]))
    sys.stderr.write('\n')

    # Set up output directory
    if not os.path.exists(outdir):
        sys.stderr.write("Creating directory: {}\n".format(outdir))
        os.mkdir(outdir)
    else:
        sys.stderr.write("{} already exists\n".format(outdir))

    # Read in GDA data tracks as pandas data frame
    # GDA data tracks file format: window_name (unique), species, chromosome, start, end, feature_values (tab separated)
    data = pd.read_csv(tracks_file, sep='\t', index_col=0)

    # Filter windows remove any N-containing windows
    if 'N_percentage' in data.columns:
        data.drop(index=data[data.N_percentage>0].index, axis=0, inplace=True)

    # Make a dataframe with only window_name index and feature values suitable for clustering
    data_to_cluster = data.drop(['species', 'start', 'end', 'chromosome'], axis=1)

    # Get UMAP embedding for the data
    sys.stderr.write("Running UMAP\n")
    embedding = run_umap(data_to_cluster, umap_n_neighbors)

    # cluster UMAP embedding using HDBSCAN
    sys.stderr.write("Running HDBSCAN\n")
    cluster_labels = run_hdbscan(embedding, leaf_size, min_samples, min_cluster_size=hdbscan_min_cluster_size)
    cluster_colours = get_cluster_cols_dict(max(set(cluster_labels)))

    # Add cluster data to main data frame
    data['cluster'] = cluster_labels

    # Write out UMAP embedding
    outfile = "umap_clustering.csv"
    write_embedding(embedding, cluster_labels, data['species'].tolist(), outfile, outdir, cluster_colours)

    # Write GFF file of clusters
    write_cluster_gff(data, outdir)

    #########################
    # Get significant features, significant feature stats per cluster, proportion of each cluster in each species and cluster means (for heatmap)
    # This step is very slow - needs optimising!
    sys.stderr.write("Determining enriched features\n")
    (sig_features, cluster_results, cluster_species_genome_prop, cluster_means) = get_significant_features(data, pvalue_cutoff)

    # Write sig features for reading by Shiny
    feature_table_filename = "feattable.csv"
    write_feature_table(cluster_results, feature_table_filename, outdir)
    # Write genome proportions
    genome_prop_outfile = "genomeprops.csv"
    write_genome_prop(cluster_species_genome_prop, genome_prop_outfile, outdir)

    # Write cluster/feature heatmap data
    cluster_heatmap_filename = 'cluster_heatmap.csv'
    write_cluster_heatmap(cluster_means, sig_features, cluster_heatmap_filename, outdir)

    ####################
    # Make data for feature histograms
    feat_hist_dict = feature_histograms(data, sig_features)
    feat_hist_filename = "feathist.json"
    write_json(feat_hist_dict, feat_hist_filename, outdir)

    #####################################
    # Species separated for these analyses
    ################
    for s in sorted(sorted(set(data['species'].tolist()))):
        # Make species-specific directory
        if not os.path.exists(outdir + '/' + s):
            os.mkdir(outdir + '/' + s)
        # Write BED file
        bed_cluster_file = write_cluster_bed(data[data.species==s], cluster_colours, outdir + '/' + s)
        # Run the script that finds which cluster junctions in the BED file occur at a different rate than what is expected by chance
        cluster_junctions_command = "cluster_junctions_fisher_test.py {} {}".format(outdir + '/' + s + '/' + bed_cluster_file, outdir + '/' + s)
        t = os.system(cluster_junctions_command)
        if t != 0:
            sys.stderr.write("Error running command {}\n".format(cluster_junctions_command))
            sys.exit(1)
        # Read back BED file of data (this could be optimised!)
        bed_data = BedFile(outdir + '/' + s + '/' + bed_cluster_file)

        #########################################
        # Generate cluster histograms with window number
        window_freqs = bed_data.cluster_histograms(cluster_position_histogram_window_number)
        # Write dict of data out to file
        cluster_position_histogram_file = "clusterpos.json"
        write_json(window_freqs, cluster_position_histogram_file, outdir + '/' + s + '/')

        #####################
        #Make heatmap of chromosome cluster composition
        comp_res = bed_data.chromosome_composition()
        chr_composition_heatmap_file = "chrcompheat.csv"
        write_chr_comp_heatmap(comp_res, chr_composition_heatmap_file, outdir + '/' + s + '/')

        ######################
        # Make JSON files for Circos data
        ######################
        circos_json_dict = bed_data.circos_json()
        write_circos_json(circos_json_dict, 'circos.json', outdir + '/' + s + '/')


if __name__ == "__main__":
    ###############
    # Parse command line arguments
    ###############
    parser = argparse.ArgumentParser(description='Genome Decomposition Analysis of windowed tracks')
    parser.add_argument("-n", "--n_neighbors", help="N neighbours argument for UMAP [13]", default=13, type=int)
    parser.add_argument("-c", "--cluster_size_cutoff", help="HDBSCAN min cluster size [200]", default=200, type=int)
    parser.add_argument("-p", "--pvalue_cutoff", help="p-value cutoff for feature enrichment in clusters [1e-20]", default=1e-20, type=float)
    parser.add_argument("-w", "--cluster_position_histogram_window_number", help="Cluster position histogram window number", default=20, type=float)
    parser.add_argument("--leaf_size", help="leaf_size setting for HDBSCAN (default: 40)", default=40, type=int)
    parser.add_argument("--min_samples", help="min_samples setting for HDBSCAN (default: None)", default=None, type=int)
    parser.add_argument("-d", "--directory", help="Output dir [gda_out]", default="gda_out", type=str)
    parser.add_argument("-m", "--metrics_file_path", help="Optional: path to a clustering metrics CSV file generated by gda_parameters.py. If this argument is used, the values for n_neighbors and cluster_size_cutoff will be read from the top row of this CSV file", default="", type=str)
    parser.add_argument("tracks", help="File of windowed GDA tracks", type=str)
    args = parser.parse_args()
    main(args.n_neighbors, args.cluster_size_cutoff, args.pvalue_cutoff, args.cluster_position_histogram_window_number, args.leaf_size, args.min_samples, args.directory, args.metrics_file_path, args.tracks)



