#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <stdexcept>
#include "boost/program_options.hpp"

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::ostream;
using std::istream;
using std::ofstream;
using std::runtime_error;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;

void run_command(string& argument_string);

void split_as_string(vector<string>& tokens, string& s, string& separators);

void split_as_double(vector<double>& tokens, string& s, string& separators);

void parse_comma_separated_pair_as_doubles(pair<double,double>& p, string& s);

double log10_sum_exp(double x1, double x2);

string join(vector <string> s, char delimiter);

variables_map parse_arguments(int argc, char* argv[], options_description options);


#endif //RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
