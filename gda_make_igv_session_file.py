#!/usr/bin/env python3
# Function(s) to generate an IGV session file for viewing GDA data

import sys
import glob
import re
import os.path
import argparse

# inputs:
# Significant features
# Clusters BED file
# FA
# Genes GFF?
# Location of BEDGraph files

# outputs:
# XML file

# Logic:
# Determine which features to show e.g. any repeats/kmers?

parser = argparse.ArgumentParser(description='Generate IGV session file from GDA results')
parser.add_argument("cluster_heatmap", help="Data file for cluster heatmap - used to get significant features", type=str)
parser.add_argument("clusters_bed", help="BED file of GDA clusters", type=str)
parser.add_argument("fa_file", help="Genome fasta file", type=str)
parser.add_argument("bedgraph_files_dir", help="Directory of bedgraph files", type=str)
parser.add_argument("-g", help="Genome annotation in GFF format", type=str)
parser.add_argument("-o", help="Output file [igv_session_gda.xml]", default="igv_session_gda.xml", type=str)
args = parser.parse_args()

if args.cluster_heatmap:
    cluster_heatmap = args.cluster_heatmap
if args.clusters_bed:
    clusters_bed = args.clusters_bed
if args.fa_file:
    fa_file = args.fa_file
if args.bedgraph_files_dir:
    bedgraph_files_dir = args.bedgraph_files_dir
if args.g:
    genes_gff_file = args.g
else:
    genes_gff_file = ""
if args.o:
    out_file = args.o

print(cluster_heatmap, clusters_bed, fa_file, bedgraph_files_dir, genes_gff_file, out_file)

class XMLString():
    '''Compose an XML string'''
    def __init__(self):
        self.xml_string = ''
        return None

    def add(self, string):
        self.xml_string = self.xml_string + string + '\n'

    def get(self):
        return self.xml_string

    def write(self, outfile):
        o = open(outfile, 'w')
        o.write(self.xml_string)
        o.close()

def get_feat_name(fn):
    o = open(fn)
    fl = o.readline()
    o.close()
    match = re.findall(r'name="([^"]*)"', fl)
    return(match[0])

def get_sig_features(datafile):
    o = open(datafile)
    l = o.readline()
    o.close()
    l = l.rstrip()
    v = l.split(',')
    return v[1:]

bedgraph_names = dict()

# Get bedgraph file paths for each feature
for bg_file in glob.glob(bedgraph_files_dir + '/*bedgraph'):
    if os.path.isfile(bg_file):
        name = get_feat_name(bg_file)
        print(name, bg_file)
        bedgraph_names[name] = bg_file

xml = XMLString()

xml.add('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
xml.add('<Session genome="{}" hasGeneTrack="false" hasSequenceTrack="true" version="8">'.format(fa_file))

xml.add('  <Resources>')
xml.add('    <Resource path="{}"/>'.format(clusters_bed))
xml.add('    <Resource path="{}"/>'.format(genes_gff_file))

# Get significant features
sig_features = get_sig_features(cluster_heatmap)
print(sig_features)

for x in sig_features:
    if x in bedgraph_names:
        xml.add('     <Resource path="{}"/>'.format(bedgraph_names[x]))
    else:
        print("Could not find a bedgraph file for {}".format(x)) 

xml.add('  </Resources>')

xml.add('</Session>')
xml.write(out_file)

print("\nPlease load {} into IGV as a session file\n".format(out_file))
