#include "BedReader.hpp"
#include <iostream>
#include <stdexcept>

using std::stoi;
using std::cout;
using std::runtime_error;


BedReader::BedReader(path bed_path){
    this->bed_path = bed_path;
    this->bed_file = ifstream(bed_path);

    if (not this->bed_file.good()){
        throw runtime_error("ERROR: bed file could not be opened: " + this->bed_path.string());
    }
}


void BedReader::read_regions(vector<Region>& regions){
    regions.emplace_back();
    while (this->next_line(regions.back())){
        regions.emplace_back();
    }

    // Deal with the empty line at the end of the file
    regions.resize(regions.size()-1);
}


bool BedReader::next_line(Region& region){
    string line;
    string start_string;
    string stop_string;

    size_t element_index = 0;
    getline(this->bed_file, line);

    if (not line.empty()){
        for (auto& character: line){
            if (character == '\t'){
                element_index++;
            }
            if (element_index == BedReader::NAME){
                region.name += character;
            }
            else if (element_index == BedReader::START){
                start_string += character;
            }
            else if (element_index == BedReader::STOP){
                stop_string += character;
            }
        }

        region.start = stoi(start_string);
        region.stop = stoi(stop_string);
    }

    return not this->bed_file.eof();
}


void BedReader::subset_by_regions_name(vector<Region>& regions, set<string> names){
    vector<Region> subset;

    for (auto& region: regions){
        if (names.count(region.name) > 0){
            subset.push_back(region);
        }
    }

    regions = move(subset);
}
