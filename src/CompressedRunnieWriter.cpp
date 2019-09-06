
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include "CompressedRunnieWriter.hpp"
#include "Miscellaneous.hpp"
#include "RunnieReader.hpp"
#include <experimental/filesystem>
#include <utility>
#include <bitset>

using boost::icl::interval_map;
using boost::icl::interval;
using boost::icl::total_enricher;
using std::experimental::filesystem::create_directories;
using std::unordered_map;
using std::make_pair;
using std::getline;
using std::bitset;
using std::cerr;


CompressedRunnieWriter::CompressedRunnieWriter(path file_path, path params_path) {
    cerr << "WRITING TO: " << file_path << '\n';

    this->sequence_file_path = file_path;
    this->index_file_path = file_path.string() + ".idx";
    this->params_path = params_path;

    // Ensure that the output directory exists
    create_directories(this->sequence_file_path.parent_path());

    this->sequence_file = ofstream(this->sequence_file_path,ofstream::binary);

    if (not this->sequence_file.is_open()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }

    this->load_parameters();
}


void CompressedRunnieWriter::build_recursive_interval_tree(){
    ///
    /// Build a recursive interval tree such that each scale interval refers to its substituent shape tree, which refers
    /// to the encoding for (shape|scale)
    ///

    // This is incremented for each cluster/2D interval
    uint8_t byte_encoding = 0;

    // Convert each vector of shape intervals into an interval tree
    for (size_t i=0; i<this->shape_intervals.size(); i++){
        pair<double,double> scale_interval = this->scale_intervals[i];
        auto shape_interval_tree = interval_map<double,uint8_t,total_enricher>();

        // Build the child tree
        for (auto& shape_interval: this->shape_intervals[i]){
            auto a = interval<double>::right_open(shape_interval.first, shape_interval.second);
            shape_interval_tree.insert(make_pair(a, byte_encoding));
            byte_encoding++;
        }

        // Add the child tree (shapes) to the parent tree (scales)
        auto b = interval<double>::right_open(scale_interval.first, scale_interval.second);
        this->recursive_interval_map.insert(make_pair(b, shape_interval_tree));
    }
}


uint8_t CompressedRunnieWriter::fetch_encoding(double scale, double shape){
    uint8_t byte = this->recursive_interval_map.find(scale)->second.find(shape)->second;

    return byte;
}


void CompressedRunnieWriter::load_parameters(){
    ifstream params_file = ifstream(this->params_path);
    vector<string> tokens;
    vector<double> bounds;
    pair<double,double> interval = {-1.0,-1.0};
    vector < pair <double,double> > intervals;
    string section;
    string line;
    string line_shapes_string;
    string line_scales_string;

    string tab_separator = "\t";
    uint64_t tab_index = 0;

    while(getline(params_file,line)){
        if (line[0] == '>'){
            section = line.substr(1,line.size());                               // Will run into the end of string
        }
        else if (section == "bounds"){
            bounds = {};
            intervals = {};
            tab_index = line.find_first_of(tab_separator);

            // Read the comma separated scale interval
            line_scales_string = line.substr(0,tab_index-1);
            parse_comma_separated_pair_as_doubles(interval, line_scales_string);
            this->scale_intervals.emplace_back(interval);

            // Read the tab separated shape bounds
            line_shapes_string = line.substr(tab_index+1,line.size());          // Will run into the end of string
            split_as_double(bounds, line_shapes_string, tab_separator);

            // Convert vector of shape bounds into vector of pairs (intervals)
            for (size_t i=0; i<bounds.size()-1; i++){
               intervals.emplace_back(bounds[i], bounds[i+1]);
            }

            // Add a vector of intervals to the shape intervals
            this->shape_intervals.emplace_back(intervals);
        }
    }

    cout << "\n\n";

    this-> build_recursive_interval_tree();
}


void CompressedRunnieWriter::write_sequence_block(RunnieSequence& sequence){
    // Write the sequence to the file
    write_string_to_binary(this->sequence_file, sequence.sequence);
}


void CompressedRunnieWriter::write_encoding_block(RunnieSequence& sequence){
    uint8_t encoding = -1;

    // Write the encodings to the file
    for (size_t i=0; i<sequence.scales.size(); i++){
        encoding = this->fetch_encoding(sequence.scales[i], sequence.shapes[i]);
        write_value_to_binary(this->sequence_file, encoding);
    }
}


void CompressedRunnieWriter::write_sequence(RunnieSequence& sequence){
    CompressedRunnieIndex index;

    // Add sequence start position to index
    index.sequence_byte_index = this->sequence_file.tellp();

    // Store the length of this sequence
    index.sequence_length = sequence.sequence.size();

    // Store the name of this sequence
    index.name = sequence.name;

    this->write_sequence_block(sequence);
    this->write_encoding_block(sequence);

    // Append index object to vector
    this->indexes.push_back(index);
}


void CompressedRunnieWriter::write_index(CompressedRunnieIndex& index){
    // Where is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_byte_index);

    // How long is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_length);

    // How long is the name of the sequence
    write_value_to_binary(this->sequence_file, index.name.size());

    // What is the name
    write_string_to_binary(this->sequence_file, index.name);
}


void CompressedRunnieWriter::write_indexes(){
    // Store the current file byte index so the beginning of the INDEX table can be located later
    uint64_t indexes_start_position = this->sequence_file.tellp();

    // Iterate all the indexes, write them to the file
    for (auto& index: this->indexes){
        write_index(index);
    }

    // Store the current file byte index so the beginning of the CHANNEL table can be located later
    uint64_t channel_metadata_start_position = this->sequence_file.tellp();

    // Write channel metadata
    write_value_to_binary(this->sequence_file, this->n_channels);
    write_value_to_binary(this->sequence_file, this->channel_size_1);

    // Write the pointer to the beginning of the index table
    write_value_to_binary(this->sequence_file, indexes_start_position);

    // Write the pointer to the beginning of the channels table
    write_value_to_binary(this->sequence_file, channel_metadata_start_position);

    cout << indexes_start_position << '\n';
    cout << channel_metadata_start_position << '\n';
}
