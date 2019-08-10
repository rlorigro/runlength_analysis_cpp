#ifndef RUNLENGTH_ANALYSIS_CPP_FASTAREADER_H
#define RUNLENGTH_ANALYSIS_CPP_FASTAREADER_H
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

using std::unordered_map;
using std::pair;
using std::vector;
using std::string;
using std::ifstream;
using std::experimental::filesystem::path;


class SequenceElement{
public:
    string name;
    string sequence;
};


class FastaIndex {
public:
    uint64_t byte_index;    // Where in the fasta is the start of the sequence
    uint64_t length;        // How long is the sequence

    FastaIndex(uint64_t byte_index, uint64_t length);
    FastaIndex();

    uint64_t size();
};


class FastaReader{
public:
    /// Attributes ///
    path file_path;
    path index_path;
    char header_symbol;
    uint64_t line_index;
    bool end_of_file;
    unordered_map <string, FastaIndex > read_indexes;

    /// Methods ///
    FastaReader(path file_path);

    // Fetch header + sequence assuming read position precedes a header
    void next_element(SequenceElement& element);

    // Generate index if needed, otherwise load index
    void index();

    // Return a copy of the read indexes
    unordered_map<string,FastaIndex> get_index();

    // Set the read_indexes attribute manually
    void set_index(unordered_map<string,FastaIndex>& read_indexes);

    // Fetch a sequence from file, and generate index first if necessary
    void fetch_sequence(SequenceElement& element, string& sequence_name);
    void fetch_sequence(SequenceElement& element, string& sequence_name, uint64_t fasta_byte_index);

private:
    /// Attributes ///
    ifstream fasta_file;

    /// Methods ///
    // Assuming the read position of the ifstream is at a header, parse the header line
    void read_next_header(SequenceElement& element);

    // Assuming the read position of the ifstream is at a sequence, parse the sequence line(s),
    // which precede the next header
    void read_next_sequence(SequenceElement& element);

    // Use htslib to make an .fai
    void build_fasta_index();

    // Read fai
    void read_fasta_index();
};


#endif
