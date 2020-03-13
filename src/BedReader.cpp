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


void BedReader::read_regions(regional_interval_map& regions){
    Region region;
    bool placeholder = true;

    while (this->next_line(region)){
        // If this contig/region hasn't been added to the map yet, add it
        if (regions.count(region.name) == 0){
            interval_map<uint64_t,bool,total_enricher> empty_interval_map;
            regions.emplace(region.name, empty_interval_map);
        }

        // Construct a Boost interval for this numeric range/coords in the BED
        auto a = interval<uint64_t>::right_open(region.start, region.stop);

        // Within this contig, add the numeric interval
        regions.at(region.name).insert(make_pair(a, placeholder));

        region = {};
    }
}


void BedReader::read_regions(vector<Region>& regions){
    // Start the vector with an empty region
    regions.emplace_back();

    while (this->next_line(regions.back())){
        // Queue an empty region for the next iteration
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
