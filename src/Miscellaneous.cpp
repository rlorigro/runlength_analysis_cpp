#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "boost/program_options.hpp"
#include <boost/tokenizer.hpp>

using std::cout;
using std::cerr;
using std::ostream;
using std::istream;
using std::string;
using std::to_string;
using std::pair;
using std::vector;
using std::runtime_error;
using std::pow;
using std::max;
using std::log10;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


// Helper function
void split_as_string(vector<string>& tokens, string& s, string& separators){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (auto& token: tok) {
        tokens.push_back(token);
    }
}


// Helper function
void split_as_double(vector<double>& tokens, string& s, string& separators){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (auto& token: tok) {
        tokens.push_back(stod(token));
    }
}


void parse_comma_separated_pair_as_doubles(pair<double,double>& p, string& s) {
    size_t comma_index = s.find_first_of(',');
    string first = s.substr(0,comma_index-1);
    string second = s.substr(comma_index+1,s.size());

    p.first = stod(first);
    p.second = stod(second);
}


size_t find_nth_character(string& s, char c, size_t n){
    size_t n_found = 0;
    size_t index = -1;

    for (size_t i=0; i<s.size(); i++){
        if (s[i] == c){
            n_found += 1;

            if (n_found == n){
                index = i;
                break;
            }
        }
    }

    return index;
}


string join(vector <string> s, char delimiter){
    string joined_string;

    for (size_t i=0; i<s.size(); i++){
        if (i < s.size()-1){
            joined_string += s[i] + delimiter;
        }
        else{
            joined_string += s[i];
        }
    }

    return joined_string;
}


double log10_sum_exp(double x1, double x2){
    ///
    /// Do a safe addition in non-log space using log values. Avoids doing multiplication in log space.
    /// sum = a + log(sum(exp(x_i-a))) for all i
    ///     where a = max(x_i) for all i

    double a = max(x1, x2);
    double sum;

    sum = a + log10(pow(10, x1-a) + pow(10, x2-a));

    return sum;
}


void print_distribution(vector<double>& distribution, uint16_t width, char character){
    ///
    /// Make a text representation of a distribution. ASSUME POSITIVE 0-1 VALUES ONLY!
    ///

    uint16_t n_characters;
    string line;

    double sum = 0;
    for (auto& y: distribution){
        sum += y;
    }

    size_t i = 0;
    for (auto& y: distribution){
        if (y < 0){
            cerr << "WARNING: Miscellaneous::print_distribution() printing negative value as 0";
            y = 0;
        }
        n_characters = uint16_t((y/sum)*width);
        line = to_string(i+1) + ":\t" + to_string(y) + "\t" + string(n_characters, character);
        cout << line << "\n";
        i++;
    }
}


variables_map parse_arguments(int argc, char* argv[], options_description options){
    try {
        // Store options in a map and apply values to each corresponding variable
        variables_map vm;
        store(parse_command_line(argc, argv, options), vm);

        // If help was specified, or no arguments given, provide help
        if (vm.count("help") || argc == 1) {
            cout << options << "\n";
            throw runtime_error("ERROR: incorrect arguments");
        }

        return vm;
    }
    // Some error in parsing arguments
    catch(boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::unknown_option> >){
        cout << options << "\n";
        throw runtime_error("ERROR: incorrect arguments");
    }
}


#endif //RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
