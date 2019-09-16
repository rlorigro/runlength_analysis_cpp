
#include "CompressedRunnieReader.hpp"


void CompressedRunnieSequence::print_encoding(){
    for (auto& element: this->encoding) {
        cout << int(element) << ',';
    }
    cout << '\n';
}


CompressedRunnieReader::CompressedRunnieReader(string file_path) {
    this->sequence_file_path = file_path;

    // Open the input file.
    this->sequence_file_descriptor = ::open(file_path.c_str(), O_RDONLY);

    // Verify it is working
    if(this->sequence_file_descriptor == -1) {
        throw runtime_error("ERROR: could not read " + file_path);
    }

    // Find file size in bytes
    this->file_length = lseek(this->sequence_file_descriptor, 0, SEEK_END);

    // Initialize remaining parameters using the file footer data
    this->read_footer();

    // Read table of contents, needed for indexed reading
    this->read_indexes();
}


size_t CompressedRunnieReader::get_read_count(){
    return this->indexes.size();
}


const string& CompressedRunnieReader::get_read_name(uint64_t read_number){
    return this->indexes.at(read_number).name;
}


uint64_t CompressedRunnieReader::get_length(uint64_t read_number){
    return this->indexes.at(read_number).sequence_length;
}


const string& CompressedRunnieReader::get_file_name(){
    return this->sequence_file_path;
}


void CompressedRunnieReader::get_sequence_data(CompressedRunnieSequence& sequence, uint64_t read_number){
    off_t byte_index = off_t(this->indexes.at(read_number).sequence_byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, sequence.sequence, this->indexes.at(read_number).sequence_length, byte_index);
    pread_vector_from_binary(this->sequence_file_descriptor, sequence.encoding, this->indexes.at(read_number).sequence_length, byte_index);
}

void CompressedRunnieReader::get_sequence_data(NamedCompressedRunnieSequence& sequence, uint64_t read_number){
    sequence.name = this->indexes.at(read_number).name;
    off_t byte_index = off_t(this->indexes.at(read_number).sequence_byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, sequence.sequence, this->indexes.at(read_number).sequence_length, byte_index);
    pread_vector_from_binary(this->sequence_file_descriptor, sequence.encoding, this->indexes.at(read_number).sequence_length, byte_index);
}


void CompressedRunnieReader::read_index_entry(CompressedRunnieIndex& index_element, off_t& byte_index){
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_byte_index, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.sequence_length, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, index_element.name_length, byte_index);
    pread_string_from_binary(this->sequence_file_descriptor, index_element.name, index_element.name_length, byte_index);
}


void CompressedRunnieReader::read_indexes(){
    off_t byte_index = off_t(this->indexes_start_position);

    while (byte_index > 0 and uint64_t(byte_index) < this->channel_metadata_start_position){
        CompressedRunnieIndex index_element;
        this->read_index_entry(index_element, byte_index);
        this->indexes.emplace_back(index_element);

        // Update the mapping of read names to their places in the vector of indexes
        auto element = make_pair(index_element.name, this->indexes.size() - 1);
        auto success = this->index_map.insert(move(element)).second;
        if (not success){
            throw runtime_error("ERROR: possible duplicate read name (" + index_element.name + ") found in runnie file: " + this->sequence_file_path);
        }
    }
}


void CompressedRunnieReader::read_channel_metadata(){
    ///
    /// Read the description of the channels that accompany the nucleotide sequence. Not currently used for initializing
    /// data vectors, since CompressedRunnieSequence.encodings specifies the data type.
    ///
    off_t byte_index = off_t(this->channel_metadata_start_position);
    pread_value_from_binary(this->sequence_file_descriptor, this->n_channels, byte_index);
    this->channel_sizes.resize(n_channels);

    for (uint64_t i=0; i<this->n_channels; i++) {
        pread_value_from_binary(this->sequence_file_descriptor, this->channel_sizes.at(i), byte_index);
    }
}


void CompressedRunnieReader::read_footer(){
    off_t byte_index = off_t(this->file_length - 2*sizeof(uint64_t));
    pread_value_from_binary(this->sequence_file_descriptor, this->indexes_start_position, byte_index);
    pread_value_from_binary(this->sequence_file_descriptor, this->channel_metadata_start_position, byte_index);

    this-> read_channel_metadata();
}
