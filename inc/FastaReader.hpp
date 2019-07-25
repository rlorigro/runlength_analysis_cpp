#ifndef RUNLENGTH_ANALYSIS_CPP_FASTAREADER_H
#define RUNLENGTH_ANALYSIS_CPP_FASTAREADER_H
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

using std::map;
using std::pair;
using std::vector;
using std::string;
using std::ifstream;
using std::experimental::filesystem::path;


struct sequence_element{
    string name;
    string sequence;
};


class FastaReader{
public:
    /// Attributes ///
    path file_path;
    path index_path;
    char header_symbol;
    uint64_t line_index;
    bool end_of_file;
    map <string,uint64_t> read_indexes;
    map <string,uint32_t> read_lengths;

    /// Methods ///
    FastaReader(path file_path);

    // Fetch header + sequence assuming read position precedes a header
    void next_element(sequence_element& element);

    // Generate index if needed, otherwise load index
    void index();

    // Return a copy of the read indexes
    map <string,uint64_t> get_index();

    // Set the read_indexes attribute manually
    void set_index(map <string,uint64_t>& read_indexes);

    // Fetch a sequence from file, and generate index first if necessary
    void fetch_sequence(sequence_element& element, string& sequence_name);
    void fetch_sequence(sequence_element& element, string& sequence_name, uint64_t fasta_byte_index);

private:
    /// Attributes ///
    ifstream fasta_file;

    /// Methods ///
    // Assuming the read position of the ifstream is at a header, parse the header line
    void read_next_header(sequence_element& element);

    // Assuming the read position of the ifstream is at a sequence, parse the sequence line(s),
    // which precede the next header
    void read_next_sequence(sequence_element& element);

    // Use htslib to make an .fai
    void build_fasta_index();

    // Read fai
    void read_fasta_index();
};


#endif
