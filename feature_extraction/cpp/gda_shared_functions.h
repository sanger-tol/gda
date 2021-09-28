/*
 * gda_shared_functions.h
 *
 *  Created on: 2 Jan 2021
 *      Author: ea
 */

#include <iostream>
//#include <fstream>
#include <string>
#include <vector>

using namespace std;


#ifndef GDA_SHARED_FUNCTIONS_H_
#define GDA_SHARED_FUNCTIONS_H_


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


#endif /* GDA_SHARED_FUNCTIONS_H_ */
