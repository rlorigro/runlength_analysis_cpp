
#ifndef RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP

#include <string>
#include <vector>
#include "Pileup.hpp"
#include "AlignedSegment.hpp"

using std::vector;
using std::string;


class RunlengthSequenceElement{
public:
//    RunlengthSequenceElement(string name, string sequence, vector<uint16_t> lengths);

    static const size_t n_channels = 1;
    string name;
    string sequence;
    vector<uint16_t> lengths;

    // Fetch the data from this sequence format associated with an index and a cigar. Alternatively, just return a null vector
    void get_read_data(vector<float>& read_data, Cigar& cigar, Coordinate& coordinate, AlignedSegment& alignment);

    // Create data to represent ambiguous sequence
    void generate_default_data_vector(vector<float>& read_data);
};

#endif //RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
