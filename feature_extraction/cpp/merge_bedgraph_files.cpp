// Program for merging bedgraph files into a TSV table

// MIT License
// 
// Copyright (c) 2020-2021 Genome Research Ltd.
//
// Author: Eerik Aunin (ea10@sanger.ac.uk)
//
// This file is a part of the Genome Decomposition Analysis (GDA) pipeline.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <asm-generic/errno-base.h>
#include <glob.h>
#include <stddef.h>
#include <sys/stat.h>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

vector<string> glob_vector(const string& pattern) {
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i) {
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}


int dir_exists(const char* const path) {
	// checks if directory exists
	// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
    struct stat info;

    int statRC = stat(path, &info);
    if (statRC != 0) {
		if (errno == ENOENT)  {return 0;} // something along the path does not exist
		if (errno == ENOTDIR) {return 0;} // something in path prefix is not a dir
		return -1;
    }
    return (info.st_mode & S_IFDIR) ? 1 : 0;
}


template<typename K, typename V>
void print_map(map<K,V> const &m) {
	// prints map keys and values
	// https://www.techiedelight.com/print-keys-values-map-cpp/
    for (auto const& pair: m) {
        cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}



map<string, string> read_fasta_as_map(string fasta_path) {
	// Input: path to a FASTA file
	// Output: map where the keys are FASTA headers and values are the corresponding sequences
	ifstream file(fasta_path);
	map<string, string> fasta_data;
	string fasta_header;
	if (file.is_open()) {
		string line;
		while (getline(file, line)) {
			int line_len = line.length();
			if (line_len > 0) {
				if (line.at(0) == '>') {
					if (line_len > 1) {
						fasta_header = line.substr(1, line_len);
						fasta_data[fasta_header] = "";
					} else {
						cerr << "Header without sequence name encountered in the FASTA file " << fasta_path << endl;
						exit(EXIT_FAILURE);
					}

				} else {
					line.erase(line.find_last_not_of(" \n\r\t") + 1);
					fasta_data[fasta_header] += line;
				}
			}
		}
		file.close();
	} else {
		perror("Cannot open file");
	}
	return fasta_data;
}

tuple<vector<string>, map<string, int>, vector<string>, vector<int>> get_window_names_map(map<string, string> fasta_data, int chunk_size) {
	// Input: 1) path to assembly FASTA file, 2) window length
    // Output: map where the keys are window names (scaffold name + window start coordinate) and the values are the indices of the windows in the assembly
	int w_counter = 0;
	vector<string> window_names;
	vector<string> scaff_names;
	vector<int> scaff_start_coords;

	map<string, int> window_names_map;
	for (auto x: fasta_data) {
		string header = x.first;
		int seq_len = x.second.length();
		for (int chunk_start = 0; chunk_start < seq_len; chunk_start += chunk_size) {
			string scaff_start = header + "_" + to_string(chunk_start);
			// cout << scaff_start << "\n";
			window_names.push_back(scaff_start);
			scaff_names.push_back(header);
			scaff_start_coords.push_back(chunk_start);
			window_names_map[scaff_start] = w_counter;
			w_counter++;
		}
	}
	return make_tuple(window_names, window_names_map, scaff_names, scaff_start_coords);
}


vector<string> split(const string& str, const string& delim) {
	// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
	// function for splitting string by a delimiter
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}


string get_bedgraph_variable_name (string line, string bedgraph_file) {
	// Input: 1) bedgraph file header line, 2) bedgraph file name
	// Output: bedgraph variable name extracted from the header line of the bedgraph file
	vector<string> split_string = split(line, "type=bedGraph name=\"");
	string bedgraph_variable = "";
	if (split_string.size() == 2) {
		vector<string> split_string2 = split(split_string[1], "\" description=");
		if (split_string2.size() == 2) {
			bedgraph_variable = split_string2[0];
		} else {
			cerr << "Failed to parse the header of bedgraph file " << bedgraph_file << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		cerr << "Failed to parse the header of bedgraph file " << bedgraph_file << endl;
		exit(EXIT_FAILURE);
	}
	replace(bedgraph_variable.begin(), bedgraph_variable.end(), ' ', '_');
	return bedgraph_variable;
}


string join(const vector<string> &lst, const string &delim) {
	// https://stackoverflow.com/questions/5689003/how-to-implode-a-vector-of-strings-into-a-string-the-elegant-way
	// function for joining string vector with delimiter
    string ret;
    for(const auto &s : lst) {
        if(!ret.empty())
            ret += delim;
        ret += s;
    }
    return ret;
}


inline bool ends_with(string const & value, string const & ending) {
	// https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
	// function for checking if a string ends with another string
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

vector<string> get_bedgraph_files_list(string in_folder, bool include_repeat_family_tracks) {
	// Input: 1) path to the pipeline output folder, 2) boolean that determines whether to include bedgraph tracks of repeat familes
	// Output: list of bedgraph files that will be used to for creating the output TSV file

	int dir_exists_check_result = dir_exists(in_folder.c_str());
	if (dir_exists_check_result != 1) {
		cerr << "Folder not found: " << in_folder << endl;
		exit(EXIT_FAILURE);
	}

	vector<string> input_folders;
	input_folders.push_back(in_folder);

	if (include_repeat_family_tracks == true) {
		string complex_repeats_folder_path = in_folder + "/complex_repeats_bedgraph";
		string simple_repeats_folder_path = in_folder + "/simple_repeats_bedgraph";
		int dir_exists_check_result_complex_repeats = dir_exists(complex_repeats_folder_path.c_str());
		if (dir_exists_check_result_complex_repeats == 1) {
			input_folders.push_back(complex_repeats_folder_path);
		} else {
			cerr << "Warning: folder not found: " << complex_repeats_folder_path << endl;
		}

		int dir_exists_check_result_simple_repeats = dir_exists(simple_repeats_folder_path.c_str());
		if (dir_exists_check_result_simple_repeats == 1) {
			input_folders.push_back(simple_repeats_folder_path);
		} else {
			cerr << "Warning: folder not found: " << simple_repeats_folder_path << endl;
		}
	}

	vector<string> bedgraph_files;
	int counter = 0;

	for(auto &folder: input_folders) {
		folder += "/*";
		vector<string> input_files = glob_vector(folder);
		for(auto const& infile: input_files) {
			//cout << infile << endl;
			if (ends_with(infile, ".bedgraph")) {
				bedgraph_files.push_back(infile);
			}
		}
		if (counter == 0) {
			sort(bedgraph_files.begin(), bedgraph_files.end());
		}
		counter++;
	}

	if (bedgraph_files.size() == 0) {
		cerr << "No bedgraph files were found in the input folder: " << in_folder << endl;
		exit(EXIT_FAILURE);
	}
	return bedgraph_files;
}


int main(int argc,char* argv[]) {
	cout.precision(16);

	string assembly_fasta_path;
	string in_folder;
	int chunk_size;
	string assembly_title;
	bool include_repeat_family_tracks;

	if (argc == 6) {
		assembly_fasta_path = argv[1];
		in_folder = argv[2];
		chunk_size = stoi(argv[3]);
		assembly_title = argv[4];
		string include_repeat_family_tracks_string = argv[5];
		if (include_repeat_family_tracks_string.compare("true") == 0) {
			include_repeat_family_tracks = true;
		} else if (include_repeat_family_tracks_string == "false") {
			include_repeat_family_tracks = false;
		} else {
			cerr << "Invalid value for the include_repeat_family_tracks variable. Expected 'true' or 'false', received '" << argv[5] << "'" << endl;
		}

	} else {
		cerr << "Error: 5 arguments expected but " << argc -1 << " received" << endl;
		exit(EXIT_FAILURE);
	}


	vector<string> bedgraph_files = get_bedgraph_files_list(in_folder, include_repeat_family_tracks);
	map<string, string> fasta_data = read_fasta_as_map(assembly_fasta_path);
	tuple<vector<string>, map<string, int>, vector<string>, vector<int>> window_names_tuple = get_window_names_map(fasta_data, chunk_size);
	vector<string> window_names = get<0>(window_names_tuple);
	map<string, int> window_names_map = get<1>(window_names_tuple);
	vector<string> scaff_names = get<2>(window_names_tuple);
	vector<int> scaff_start_coords = get<3>(window_names_tuple);

	int window_names_map_size = window_names_map.size();
	int bedgraph_files_count = bedgraph_files.size();
	double bg_array[window_names_map_size][bedgraph_files_count] = {{0}};
	vector<string> bedgraph_variables;

	int file_counter = 0;
	for (auto &bedgraph_file: bedgraph_files) {

		string bedgraph_variable = "";
		ifstream file(bedgraph_file);

		if (file.is_open()) {
		    string line;

		    int line_counter = 0;
		    while (getline(file, line)) {
		    	if (line_counter == 0) {
		    		bedgraph_variable = get_bedgraph_variable_name(line, bedgraph_file);
					bedgraph_variables.push_back(bedgraph_variable);
		    	} else {

		    		vector<string> split_line = split(line, " ");

		    		if (split_line.size() != 4) {
		    			cerr << "Failed to parse the bedgraph file " << bedgraph_file << " at line " << line_counter << ":" << endl;
		    			cerr << line << endl;
		    			exit(EXIT_FAILURE);
		    		} else {

		    			string window = split_line[0] + "_" + split_line[1];
		    			double window_value = strtod(split_line[3].c_str(), NULL);

		    			if (window_names_map.find(window) == window_names_map.end()) {
							cerr << "Window " << window << " from bedgraph file " << bedgraph_file << " not found in window_names_map derived from the input FASTA file" << endl;
                            exit(EXIT_FAILURE);
		    			} else {
		    				int window_array_pos = window_names_map[window];
		    				bg_array[window_array_pos][file_counter] = window_value;
		    			}
		    		}
		    	}
		    	line_counter++;
		    }
		    file.close();
		}
		file_counter++;

	}

	string header_row = "window\tstart\tend\tspecies\tchromosome\t" + join(bedgraph_variables, "\t");
	cout << header_row << endl;

	for (int i = 0; i < window_names_map_size; i++) {
		string window_name = window_names[i];
		string scaff = scaff_names[i];

		int start_coord = scaff_start_coords[i];
		int end_coord = start_coord + chunk_size - 1;

		cout << window_name << "\t" << start_coord << "\t" << end_coord << "\t" << assembly_title << "\t" << scaff << "\t";

		for (int j = 0; j < bedgraph_files_count; j++) {
			cout << bg_array[i][j];
			if (j < bedgraph_files_count -1) cout << "\t";
		}
		cout << endl;
	}
	return 0;
}

