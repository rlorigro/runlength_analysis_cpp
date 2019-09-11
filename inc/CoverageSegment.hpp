
#ifndef RUNLENGTH_ANALYSIS_COVERAGESEGMENT_HPP
#define RUNLENGTH_ANALYSIS_COVERAGESEGMENT_HPP


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

    /// Methods ///
    void print();
};

#endif //RUNLENGTH_ANALYSIS_COVERAGEELEMENT_HPP
