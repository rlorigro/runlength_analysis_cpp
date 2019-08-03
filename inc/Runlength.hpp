#ifndef RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
#define RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP

#include "FastaReader.hpp"
#include <vector>

using std::vector;


class RunlengthSequenceElement{
public:
    string name;
    string sequence;
    vector<uint16_t> lengths;
};

void runlength_encode(RunlengthSequenceElement& runlength_sequence, SequenceElement& sequence);


#endif //RUNLENGTH_ANALYSIS_CPP_RUNLENGTH_HPP
