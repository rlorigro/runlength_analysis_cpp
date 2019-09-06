#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include "boost/program_options.hpp"

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::ostream;
using std::istream;
using std::ofstream;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;

void run_command(string& argument_string);

void split_as_string(vector<string>& tokens, string& s, string& separators);

void split_as_double(vector<double>& tokens, string& s, string& separators);

void parse_comma_separated_pair_as_doubles(pair<double,double>& p, string& s);

string join(vector <string> s, char delimiter);

variables_map parse_arguments(int argc, char* argv[], options_description options);


void write_string_to_binary(ostream& s, string& stringaling);

void read_string_from_binary(istream& s, string& v, uint64_t length);


template<class T> void write_vector_to_binary(ostream& s, const vector<T>& v){
    ///
    /// Without worrying about size conversions, write any vector to a file using ostream.write
    ///

    s.write(reinterpret_cast<const char*>(v.data()), v.size()*sizeof(T));
}


template<class T> void write_value_to_binary(ostream& s, T v){
    ///
    /// Without worrying about size conversions, write any value to a file using ostream.write
    ///

    auto v_temp = v;
    s.write(reinterpret_cast<const char*>(&v_temp), sizeof(T));
}


template<class T> void read_value_from_binary(istream& s, T& v){
    ///
    /// Without worrying about size conversions, write any value to a file using ostream.write
    ///
    cout << "Reading value size of: " << sizeof(T) << " at position: " << s.tellg() << '\n';
    s.read(reinterpret_cast<char*>(&v), sizeof(T));
}


#endif //RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
