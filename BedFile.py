#!/usr/bin/env python3
"""
Class for representing BED files
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

class BedFile():
    ''' Parse bed file to get (amongst other things) composition of chromosomes in terms of features'''
    def __init__(self, bed_file):
        self.chromosomes = set()
        self.all_features = set()
        self.features = dict() # dict of chromosome of start positions of list (end, feature_type, color)
        self.chr_lengths = dict()
        self.features_per_chromosome = dict() # Count of the total number of features per chromosome to allow exclusion of chromosomes with too few features for windowing

        b = open(bed_file)
        for x in b.readlines():
        #tarseq_0_pilon  0       5000    set_5   0       +       0       5000    #E76BF3
            x = x.rstrip()
            if x.startswith('#'):
                continue
            v = x.split('\t')
            self.chromosomes.add(v[0])
            self.all_features.add(v[3])
            if v[0] not in self.features:
                self.features[v[0]] = dict()
            if v[1] not in self.features[v[0]]:
                self.features[v[0]][v[1]] = list()
            self.features[v[0]][v[1]] = [v[2], v[3], v[8]]

            # Record chromosome length i.e. rightmost feature end (which is not guaranteed to be the chromosome end!)
            if v[0] not in self.chr_lengths:
                self.chr_lengths[v[0]] = int(v[2])
            elif self.chr_lengths[v[0]] < int(v[2]):
                self.chr_lengths[v[0]] = int(v[2])

            # Record number of features per chromosome
            if v[0] not in self.features_per_chromosome:
                self.features_per_chromosome[v[0]] = 0
            self.features_per_chromosome[v[0]] = self.features_per_chromosome[v[0]] + 1

        b.close()

    def circos_json(self):
        json_dict = dict()
        json_dict['genome'] = list()
        for c in self.chr_lengths:
            json_dict['genome'].append({'id': c, 'label': c, 'color': 'red', 'len': self.chr_lengths[c]})
        json_dict['clusters'] = list()
        for c in self.features:
            for s in self.features[c]:
                json_dict['clusters'].append({'name': c + str(s), 'block_id': c, 'start': s, 'end': self.features[c][s][0], 'cluster': self.features[c][s][1], 'color': self.features[c][s][2]})
        return(json_dict)

    def chromosome_composition(self):
        sum = dict()
        total = dict()
        for c in self.features:
            for s in self.features[c]:
                if c not in sum:
                    sum[c] = dict()
                if self.features[c][s][1] not in sum[c]:
                    sum[c][self.features[c][s][1]] = 0
                sum[c][self.features[c][s][1]] = sum[c][self.features[c][s][1]] + (int(self.features[c][s][0]) - int(s) + 1)
                if c not in total:
                    total[c] = 0
                total[c] = total[c] + (int(self.features[c][s][0]) - int(s) + 1)

        comp_res = dict()

        for c in sum:
            for f in sum[c]:
                if c not in comp_res:
                    comp_res[c] = dict()
                freq = (sum[c][f] / total[c])
                comp_res[c][f] = freq

        return comp_res

    def cluster_histograms(self, windows = 100):
        ''' Generate data describing frequency of each cluster over each chromosome in N windows
        This will be used for comparing patterns of each cluster over different chromosomes '''

        # n.b. here feature means cluster feature i.e. -1, 0, 1, 2, etc.
        feat_patterns = dict() # Dict of clusters, of dict of chromosomes of list of values per window

        for c in self.features:
            # Exclude chromosomes with fewer features than we are looking at windows
            if c in self.features_per_chromosome and self.features_per_chromosome[c] < windows:
                continue

            win_length = self.chr_lengths[c] / windows
            win_end = win_length
            win_num = 0
            win_counts = dict()
            total_win_feat = 0


            for s in sorted([int(x) for x in self.features[c]]):
                s = str(s)
                end = self.features[c][s][0]
                feat = self.features[c][s][1]
                if int(s) > win_end:

                    for f in self.all_features:
                        if f not in win_counts:
                            win_counts[f] = 0
                        freq = 0
                        if total_win_feat > 0:
                            freq = win_counts[f] / total_win_feat
                        if f not in feat_patterns:
                            feat_patterns[f] = dict()
                        if c not in feat_patterns[f]:
                            feat_patterns[f][c] = list()
                        feat_patterns[f][c].append(freq)

                    win_num = win_num + 1
                    win_end = win_end + win_length
                    win_counts = dict()
                    total_win_feat = 0
                else:
                    if feat not in win_counts:
                        win_counts[feat] = 0
                    win_counts[feat] = win_counts[feat] + 1
                    total_win_feat = total_win_feat + 1
            # Clean up final window
            for f in self.all_features:
                if f not in win_counts:
                    win_counts[f] = 0
                freq = 0
                if total_win_feat > 0:
                    freq = win_counts[f] / total_win_feat
                if f not in feat_patterns:
                    feat_patterns[f] = dict()
                if c not in feat_patterns[f]:
                    feat_patterns[f][c] = list()
                feat_patterns[f][c].append(freq)

        return feat_patterns

