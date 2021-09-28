#!/usr/bin/env python3 
"""
Script for converting Tandem Repeats Finder repeat density values to bedgraph format
Argument1: path to an output file of trf_repeat_density_sliding_window.py
Output (STDOUT): input file converted to bedgraph format
"""
import sys
import os.path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf


print('track type=bedGraph name="tandem_repeats_fraction" description="tandem_repeats_fraction" visibility=full color=0,0,255 altColor=0,100,200 priority=20')

in_path = sys.argv[1]
in_data = gpf.l(in_path)
in_data = in_data[1: len(in_data)]
for line in in_data:
    split_line = line.split()
    out_line = " ".join([split_line[0], split_line[1], str(int(split_line[2]) + 1), split_line[3]])
    print(out_line)
    

