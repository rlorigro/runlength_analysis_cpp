#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include "boost/program_options.hpp"
#include <boost/tokenizer.hpp>

using std::cout;
using std::string;
using std::vector;
using std::runtime_error;
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

    for (string token: tok) {
        tokens.push_back(token);
    }
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
