#include "MarginPolishReader.hpp"
#include <experimental/filesystem>
#include <iostream>
#include <algorithm>

using std::cout;
using std::replace;
using std::experimental::filesystem::directory_iterator;
using std::experimental::filesystem::path;


MarginPolishReader::MarginPolishReader(path directory_path){
    this->directory_path = directory_path;
}


void MarginPolishReader::index(){
    path filename_prefix;
    string read_name;

    for (const path& file_path: directory_iterator(this->directory_path)){
        filename_prefix = file_path.filename();
        filename_prefix.replace_extension("");
        read_name = filename_prefix;
        replace(read_name.begin(), read_name.end(), '.', '_');

        this->read_paths[read_name] = file_path;
    }
}


map<string,path> MarginPolishReader::get_index(){
    return this->read_paths;
}

