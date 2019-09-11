
#ifndef RUNLENGTH_ANALYSIS_COVERAGEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_COVERAGEELEMENT_HPP

#include "Base.hpp"
#include <string>
#include <vector>
#include <iostream>


class CoverageElement{
public:
    /// Attributes ///
    string base;
    uint16_t length;
    bool reversal;
    double weight;

    CoverageElement(string base, uint16_t length, bool reversal, double weight);
    static const string reversal_string;

    /// Methods ///
    string to_string();
    uint8_t get_base_index();
    bool is_conventional_base();
};

#endif //RUNLENGTH_ANALYSIS_COVERAGEELEMENT_HPP
