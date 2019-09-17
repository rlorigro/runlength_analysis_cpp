#ifndef RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
#define RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP

#include "FastaReader.hpp"
#include "Matrix.hpp"
#include <vector>

using std::vector;


class RunlengthSequenceElement{
public:
//    RunlengthSequenceElement(string name, string sequence, vector<uint16_t> lengths);

    string name;
    string sequence;
    vector<uint16_t> lengths;
};


template<class T> void runlength_encode(RunlengthSequenceElement& runlength_sequence, T& sequence);


path runlength_encode_fasta_file(path input_file_path,
                                 unordered_map <string,RunlengthSequenceElement>& runlength_sequences,
                                 path output_dir,
                                 bool store_in_memory,
                                 uint16_t max_threads);


//template <typename T> void measure_runlength_distribution_from_coverage_data(path input_directory,
//                                                       path reference_fasta_path,
//                                                       path output_directory,
//                                                       uint16_t max_runlength,
//                                                       uint16_t max_threads);

void measure_runlength_distribution_from_marginpolish(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);

void measure_runlength_distribution_from_shasta(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);


void measure_runlength_distribution_from_fasta(path reads_fasta_path,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);


void measure_runlength_distribution_from_runnie(path runnie_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_runlength,
        uint16_t max_threads);


void get_vector_from_index_map(vector< pair <string,FastaIndex> >& items, unordered_map<string,FastaIndex>& map_object);



/// TEMPLATE BASED METHODS ///

template<class T> void runlength_encode(RunlengthSequenceElement& runlength_sequence, T& sequence){
    runlength_sequence = {};

    // First just copy the name
    runlength_sequence.name = sequence.name;

    char current_character = 0;

    // Iterate each run and append/update base/length for their corresponding vectors
    for (auto& character: sequence.sequence){
        if (tolower(character) != tolower(current_character)){
            runlength_sequence.sequence += character;
            runlength_sequence.lengths.push_back(1);
        }
        else{
            runlength_sequence.lengths.back()++;
        }

        current_character = character;
    }
}


#endif //RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
