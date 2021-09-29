// Program for downsampling the TSV file that has been produced from merging bedgraph files

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <gda_shared_functions.h>

using namespace std;


void process_scaff_buffer(vector<string> scaff_buffer) {
	// Input: a vector that contains rows from the input table that will form one row in the output table
	// Output: output table row that has been generated by averaging input row

	vector<string> split_line_0 = split(scaff_buffer.at(0), "\t");
	string window_name = split_line_0.at(0);
	string species = split_line_0.at(3);
	int chunk_size = stoi(split_line_0.at(2)) - stoi(split_line_0.at(1)) + 1;
	int window_start = stoi(split_line_0.at(1));
	int scaff_buffer_len = scaff_buffer.size();
	int window_end = window_start + (chunk_size * scaff_buffer_len) - 1;
	string current_scaff_name = split_line_0.at(4);
	int ncol = split_line_0.size() - 5;
	double buffer_mat[scaff_buffer_len][ncol] = {{0}};
	int row_counter = 0;
	for (auto &item: scaff_buffer) {
		vector<string> split_item = split(item, "\t");
		split_item = std::vector<string>(split_item.begin() + 5, split_item.end());
		int col_counter = 0;
		for (auto &row_element: split_item) {
			double row_element_float = 0;
			try {
				row_element_float = stod(row_element);
			} catch (const invalid_argument&) {
				cerr << "Failed to convert string to double: " << row_element << endl;
				cerr << "Argument is invalid" << endl;
			} catch (const out_of_range&) {
				cerr << "Argument is out of range for a double: " << row_element << endl;
				throw;
			}
			buffer_mat[row_counter][col_counter] = row_element_float;
			col_counter++;
		}
		row_counter++;
	}

	cout << window_name << "\t" << window_start << "\t" << window_end << "\t" << species << "\t" << current_scaff_name << "\t";
	for(int col_index = 0; col_index < ncol; col_index++) {
		double col_sum = 0;
		for(int row_index = 0; row_index < scaff_buffer_len; row_index++) {
			col_sum += buffer_mat[row_index][col_index];
		}
		double col_mean = col_sum/scaff_buffer_len;
			cout << col_mean;
		if (col_index < ncol - 1) {
			cout << "\t";
		}
	}
	cout << endl;

}


int main(int argc,char* argv[]) {
	cout.precision(16);
	string in_path;
	int downsampling_factor;

	if (argc == 3) {
		in_path = argv[1];
		downsampling_factor = stoi(argv[2]);

	} else {
		cerr << "Error: 2 arguments expected but " << argc -1 << " received" << endl;
		exit(EXIT_FAILURE);
	}

	int max_scaff_lines_counter = downsampling_factor - 1;

	int scaff_lines_counter = 0;
	string previous_scaff_name = "";
	vector<string> scaff_buffer;

	ifstream file(in_path);
	if (file.is_open()) {
		string line;

		int counter = 0;
		while (getline(file, line)) {
			if (counter == 0) {
				cout << line << endl;
			} else if (counter > 0) {
				vector<string> split_line = split(line, "\t");
				string scaff_name = split_line.at(4);
				if (previous_scaff_name.compare("") == 0) {
					previous_scaff_name = scaff_name;
				}
				if (previous_scaff_name.compare(scaff_name) == 0) {
					if (scaff_lines_counter == max_scaff_lines_counter) {
						scaff_buffer.push_back(line);
						process_scaff_buffer(scaff_buffer);
						scaff_buffer.clear();
						scaff_lines_counter = 0;
					} else {
						scaff_buffer.push_back(line);
						scaff_lines_counter++;
					}

				} else {
					if (scaff_buffer.size() > 0) {
						process_scaff_buffer(scaff_buffer);
					}
					scaff_buffer.clear();
					scaff_buffer.push_back(line);
					scaff_lines_counter = 1;
					previous_scaff_name = scaff_name;
				}
			}
			counter++;
		}
		if (scaff_buffer.size() > 0) {
			process_scaff_buffer(scaff_buffer);
		}
		file.close();
	} else {
		perror("Cannot open file");
	}
}


