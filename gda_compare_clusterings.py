#!/usr/bin/env python3
"""
Compare two different GDA clusterings of the same data by displaying them as a graph/network
"""
# MIT License
# 
# Copyright (c) 2020-2021 Genome Research Ltd.
#
# Author: Adam Reid (ar11@sanger.ac.uk)
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

import sys
import networkx as nx
import matplotlib.pyplot as plt

c1_file = sys.argv[1]
c2_file = sys.argv[2]

min_shared = 0
edge_weight_factor = 10 # scale the edge weights

def get_clusters(file):
    '''Given a gff file 'clusters.gff' from GDA, extract the windows and their cluster identity'''
    f = open(file)
    results = dict()
    for x in f.readlines():
        x = x.rstrip()
        v = x.split('\t')
        name = v[0] + '_' + v[3]
        id, other = v[8].split(';')
        id = id.replace('ID=', '')
        results[name] = id
    return results

c1 = get_clusters(c1_file)
c2 = get_clusters(c2_file)

# Check keys are the same
if c1.keys() != c2.keys():
    print("Clusterings have different sets of windows - are they from the same genome with the same window size?")
    exit()

# Collect data
edge_weights = dict()

Anodes = set()
Bnodes = set()

for n in c1:
    Anode = 'A' + c1[n]
    Anodes.add(Anode)
    Bnode = 'B' + c2[n]
    Bnodes.add(Bnode)

    # Set up dictionary
    if Anode not in edge_weights:
        edge_weights[Anode] = dict()
    if Bnode not in edge_weights[Anode]:
        edge_weights[Anode][Bnode] = 0

    #upgrade weight
    edge_weights[Anode][Bnode] = edge_weights[Anode][Bnode] + 1

# Generate network
G = nx.Graph()
i = 0
# Sorted and reverse list of nodes
for a in sorted(Anodes)[::-1]:
    G.add_node(a,pos=(1,i))
    print (a, 1, i)
    i = i + 1

j = 0
# Sorted and reverse list of nodes
for b in sorted(Bnodes)[::-1]:
    G.add_node(b,pos=(3,j))
    print (b, 3, j)
    j = j + 1

min_edge_weight = 100000000
max_edge_weight = 0

for a in edge_weights:
    for b in edge_weights[a]:
        if edge_weights[a][b] < min_shared:
            continue
        G.add_edge(a, b, weight=edge_weights[a][b])
        if min_edge_weight > edge_weights[a][b]:
            min_edge_weight = edge_weights[a][b]
        if max_edge_weight < edge_weights[a][b]:
            max_edge_weight = edge_weights[a][b]

pos=nx.get_node_attributes(G,'pos')	# Simple layout

nx.draw_networkx_nodes(G, pos, node_size=700)

for (u, v, d) in G.edges(data=True):
    norm_width = ((d['weight'] - min_edge_weight) * edge_weight_factor) / max_edge_weight
    print(u, v, d, norm_width)
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width = norm_width)

nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")

plt.axis("off")

plt.show()


