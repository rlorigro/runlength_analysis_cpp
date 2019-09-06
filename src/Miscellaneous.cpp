#ifndef RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
#define RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>
#include "boost/program_options.hpp"
#include <boost/tokenizer.hpp>

using std::cout;
using std::ostream;
using std::istream;
using std::string;
using std::pair;
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


void write_string_to_binary(ostream& s, string& stringaling){
    ///
    /// Without worrying about size conversions, write any string to a file using ostream.write
    ///

    s.write(reinterpret_cast<const char*>(stringaling.data()), stringaling.size());
}


void read_string_from_binary(istream& s, string& stringaling, uint64_t length){
    ///
    /// Without worrying about size conversions, read any value to a file using ostream.write
    ///
    cout << "Reading value size of: " << length << " at position: " << s.tellg() << '\n';

    stringaling.resize(length);
    s.read(reinterpret_cast<char*>(stringaling.data()), length);
}


void pread_bytes(int file_descriptor, char* buffer_pointer, size_t bytes_to_read, off_t& byte_index){
    while (bytes_to_read) {
        const ssize_t byte_count = ::pread(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
        if (byte_count <= 0) {
            throw runtime_error("Error " + std::to_string(errno) + " while reading: " + string(::strerror(errno)));
        }
        bytes_to_read -= byte_count;
        buffer_pointer += byte_count;
        byte_index += byte_count;
    }
}

void pread_string_from_binary(int file_descriptor, string& s, uint64_t length, off_t& byte_index){
    ///
    /// Same as the non-p version of this function, but instead is implemented with Linux pread, which is threadsafe
    ///

    s.resize(length);

    size_t bytes_to_read = length;
    char* buffer_pointer = reinterpret_cast<char*>(s.data());

    pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


#endif //RUNLENGTH_ANALYSIS_CPP_MISCELLANEOUS_H
