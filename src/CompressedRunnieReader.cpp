
#include "CompressedRunnieReader.hpp"
#include "CompressedRunnieWriter.hpp"
#include "Miscellaneous.hpp"
#include <fcntl.h>


void CompressedRunnieSequence::print_encoding(){
    for (auto& element: this->encoding) {
        cout << int(element) << ',';
    }
    cout << '\n';
}


CompressedRunnieReader::CompressedRunnieReader(path file_path) {
    this->sequence_file_path = file_path;
//    this->params_path = params_path;

    // Open the input file.
    this->sequence_file_descriptor = ::open(file_path.c_str(), O_RDONLY);

    // Verify it is working
    if(this->sequence_file_descriptor == -1) {
        throw runtime_error("ERROR: could not read " + file_path.string());
    }

    // Find file size in bytes
    this->file_length = lseek(this->sequence_file_descriptor, 0, SEEK_END);
}


void CompressedRunnieReader::read_sequence(CompressedRunnieSequence& sequence, CompressedRunnieIndex& index_element){
    off_t byte_index = index_element.sequence_byte_index;
    pread_string_from_binary(this->sequence_file_descriptor, sequence.sequence, index_element.sequence_length, byte_index);
    pread_vector_from_binary(this->sequence_file_descriptor, sequence.encoding, index_element.sequence_length, byte_index);
}


void CompressedRunnieReader::read_index_entry(CompressedRunnieIndex& index_element, off_t& byte_index){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_byte_index, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.name_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.name, index_element.name_length, byte_index);
}


void CompressedRunnieReader::read_indexes(){
    off_t byte_index = this->indexes_start_position;

    while (byte_index > 0 and uint64_t(byte_index) < this->channel_metadata_start_position){
        CompressedRunnieIndex index_element;
        this->read_index_entry(index_element, byte_index);
        this->indexes.emplace_back(index_element);

        // Update the mapping of read names to their places in the vector of indexes
        auto success = this->index_map.try_emplace(index_element.name, this->indexes.size() - 1).second;
        if (not success){
            throw runtime_error("ERROR: possible duplicate read name (" + index_element.name + ") found in runnie file: " + this->sequence_file_path.string());
        }
    }
}


void CompressedRunnieReader::read_footer(){
    off_t byte_index = this->file_length - 2*sizeof(uint64_t);
    pread_value_from_binary(this->sequence_file_descriptor, this->indexes_start_position, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, this->channel_metadata_start_position, byte_index);
}
