
#ifndef RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP

#include "Pileup.hpp"
#include "AlignedSegment.hpp"
#include <vector>
#include <string>

using std::vector;
using std::string;


class SequenceElement{
public:
    string name;
    string sequence;
    static const size_t n_channels = 1;

    // Fetch the data from this sequence format associated with an index and a cigar. Alternatively, just return a null vector
    void get_read_data(vector<float>& read_data, Cigar& cigar, Coordinate& coordinate);

    // Create data to represent ambiguous sequence
    void generate_default_data_vector(vector<float>& read_data);

//    SequenceElement();
};


#endif //RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP
