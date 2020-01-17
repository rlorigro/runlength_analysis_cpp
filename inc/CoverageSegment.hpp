
#ifndef RUNLENGTH_ANALYSIS_COVERAGESEGMENT_HPP
#define RUNLENGTH_ANALYSIS_COVERAGESEGMENT_HPP


#include "RunlengthSequenceElement.hpp"
#include "CoverageElement.hpp"
#include "Base.hpp"
#include <string>
#include <vector>
#include <iostream>

using std::vector;
using std::string;
using std::cout;


class CoverageSegment{
public:
    /// Attributes ///
    string name;
    vector <vector <CoverageElement> > coverage_data;
    string sequence;
    vector<uint16_t> lengths;
    vector<uint16_t> n_coverage;
    vector<bool> is_vertex;

    /// Methods ///
    void print();
};

#endif //RUNLENGTH_ANALYSIS_COVERAGEELEMENT_HPP
