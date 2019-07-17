#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <vector>
#include <string>
#include "boost/program_options.hpp"

using std::string;
using std::vector;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;

string join(vector <string> s, char delimiter);

variables_map parse_arguments(int argc, char* argv[], options_description options);

#endif //RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
