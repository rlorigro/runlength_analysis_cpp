
#ifndef RUNLENGTH_ANALYSIS_RUNNIESEQUENCEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_RUNNIESEQUENCEELEMENT_HPP

#include <string>
#include <vector>
#include "AlignedSegment.hpp"

using std::vector;
using std::string;


class RunnieSequenceElement{
public:
    /// Attributes ///
    size_t n_channels = 3;

    string name;
    string sequence;
    vector <float> scales;
    vector <float> shapes;

    /// Methods ///
    RunnieSequenceElement()=default;

    // Fetch the data from this sequence format associated with an index and a cigar. Alternatively, just return a null vector
    void get_read_data(vector<float>& read_data, Cigar& cigar, Coordinate& coordinate);

    // Create data to represent ambiguous sequence
    void generate_default_data_vector(vector<float>& read_data);

};


#endif //RUNLENGTH_ANALYSIS_RUNNIESEQUENCEELEMENT_HPP
