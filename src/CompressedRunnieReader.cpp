
#include "CompressedRunnieReader.hpp"
#include "CompressedRunnieWriter.hpp"
#include "Miscellaneous.hpp"


CompressedRunnieReader::CompressedRunnieReader(path file_path) {
    this->sequence_file_path = file_path;
//    this->params_path = params_path;

    this->sequence_file = ifstream(this->sequence_file_path);

    // Find file size in bytes
    this->sequence_file.seekg(0, this->sequence_file.end);
    this->file_length = this->sequence_file.tellg();

    if (not this->sequence_file.good()) {
        throw runtime_error("ERROR: could not write file " + file_path.string());
    }
}


void CompressedRunnieReader::read_index_entry(CompressedRunnieIndex& index_element){
    read_value_from_binary(this->sequence_file, index_element.sequence_byte_index);
    read_value_from_binary(this->sequence_file, index_element.sequence_length);
    read_value_from_binary(this->sequence_file, index_element.name_length);
    read_string_from_binary(this->sequence_file, index_element.name, index_element.name_length);
    cout << "INDEX DATA: " << index_element.sequence_byte_index << " " << index_element.sequence_length << " " << index_element.name << '\n';
}


void CompressedRunnieReader::read_footer(){
    this->sequence_file.seekg(-2*sizeof(uint64_t), this->sequence_file.end);
    read_value_from_binary(this->sequence_file, this->indexes_start_position);

    this->sequence_file.seekg(-sizeof(uint64_t), this->sequence_file.end);
    read_value_from_binary(this->sequence_file, this->channel_metadata_start_position);

    cout << "channel_metadata_start_position: " << this->channel_metadata_start_position << '\n';
    cout << "indexes_start_position: " << this->indexes_start_position << '\n';
}


void CompressedRunnieReader::read_indexes(){
    this->sequence_file.seekg(this->indexes_start_position);

    while (this->sequence_file.tellg() > 0 and uint64_t(this->sequence_file.tellg()) < this->channel_metadata_start_position){
        cout << this->sequence_file.tellg() << '\n';

        CompressedRunnieIndex index;
        this->read_index_entry(index);
        this->indexes.emplace_back(index);

        // Update the mapping of read names to their places in the vector of indexes
        auto success = this->index_map.try_emplace(index.name, this->indexes.size()-1).second;
        if (not success){
            throw runtime_error("ERROR: possible duplicate read name (" + index.name + ") found in runnie file: " + this->sequence_file_path.string());
        }
    }
}
