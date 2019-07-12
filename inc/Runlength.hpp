#ifndef RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
#define RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP

#include "FastaReader.hpp"
#include <vector>

using std::vector;


struct runlength_sequence_element{
    string name;
    string sequence;
    vector<uint16_t> lengths;
};

void runlength_encode(runlength_sequence_element& runlength_sequence, sequence_element& sequence);


#endif //RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
