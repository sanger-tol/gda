#!/bin/bash
cpp_files_path=$1
mkdir -p ${cpp_files_path}/bin
g++ -std=c++11 ${cpp_files_path}/downsample_merged_bedgraph_file.cpp -I ${cpp_files_path}/ -o ${cpp_files_path}/bin/downsample_merged_bedgraph_file
g++ -std=c++11 ${cpp_files_path}/merge_bedgraph_files.cpp -I ${cpp_files_path}/ -o ${cpp_files_path}/bin/merge_bedgraph_files
chmod +x ${cpp_files_path}/bin/*

