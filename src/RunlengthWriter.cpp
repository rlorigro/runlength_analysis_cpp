
#include "RunlengthWriter.hpp"
#include "BinaryIO.hpp"
#include <experimental/filesystem>
#include <bitset>

using std::experimental::filesystem::create_directories;
using std::unordered_map;
using std::make_pair;
using std::getline;
using std::bitset;
using std::cerr;

const vector<uint64_t> RunlengthWriter::channel_sizes = {sizeof(uint16_t)};

ostream& operator<<(ostream& s, RunlengthIndex& index) {
    s << "sequence_length:\t" << index.sequence_length << '\n';
    s << "sequence_byte_index:\t" << index.sequence_byte_index << '\n';
    s << "name_length:\t" << index.name_length << '\n';
    s << "name:\t" << index.name;

    return s;
}


RunlengthWriter::RunlengthWriter(path file_path) {
    this->sequence_file_path = file_path;

    // Ensure that the output directory exists
    create_directories(this->sequence_file_path.parent_path());

    this->sequence_file = ofstream(this->sequence_file_path,ofstream::binary);

    if (not this->sequence_file.is_open()){
        throw runtime_error("ERROR: could not open file " + file_path.string());
    }
}


void RunlengthWriter::write_sequence_block(RunlengthSequenceElement& sequence){
    // Write the sequence to the file
    write_string_to_binary(this->sequence_file, sequence.sequence);
}


void RunlengthWriter::write_length_block(RunlengthSequenceElement& sequence){
    // Write the encodings to the file
    for (auto& length: sequence.lengths){
        write_value_to_binary(this->sequence_file, length);
    }
}


void RunlengthWriter::write_sequence(RunlengthSequenceElement& sequence){
    if (sequence.sequence.empty()){
        throw runtime_error("ERROR: empty sequence provided to RunlengthWriter: " + sequence.name);
    }

    RunlengthIndex index;

    // Add sequence start position to index
    index.sequence_byte_index = this->sequence_file.tellp();

    // Store the length of this sequence
    index.sequence_length = sequence.sequence.size();

    // Store the name of this sequence
    index.name = sequence.name;

    this->write_sequence_block(sequence);
    this->write_length_block(sequence);

    // Append index object to vector
    this->indexes.push_back(index);
}


void RunlengthWriter::write_index(RunlengthIndex& index){
    // Where is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_byte_index);

    // How long is the sequence
    write_value_to_binary(this->sequence_file, index.sequence_length);

    // How long is the name of the sequence
    write_value_to_binary(this->sequence_file, index.name.size());

    // What is the name
    write_string_to_binary(this->sequence_file, index.name);
}


void RunlengthWriter::write_indexes(){
    // Store the current file byte index so the beginning of the INDEX table can be located later
    uint64_t indexes_start_position = this->sequence_file.tellp();

    // Iterate all the indexes, write them to the file
    for (auto& index: this->indexes){
        write_index(index);
    }

    // Store the current file byte index so the beginning of the CHANNEL table can be located later
    uint64_t channel_metadata_start_position = this->sequence_file.tellp();

    // Write channel metadata
    write_value_to_binary(this->sequence_file, RunlengthWriter::n_channels);
    write_value_to_binary(this->sequence_file, RunlengthWriter::channel_sizes[RunlengthWriter::LENGTH]);

    // Write the pointer to the beginning of the index table
    write_value_to_binary(this->sequence_file, indexes_start_position);

    // Write the pointer to the beginning of the channels table
    write_value_to_binary(this->sequence_file, channel_metadata_start_position);
}
