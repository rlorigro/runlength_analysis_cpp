
#ifndef RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP

#include <string>
#include <vector>

using std::vector;
using std::string;


class RunlengthSequenceElement{
public:
//    RunlengthSequenceElement(string name, string sequence, vector<uint16_t> lengths);

    string name;
    string sequence;
    vector<uint16_t> lengths;
};

#endif //RUNLENGTH_ANALYSIS_RUNLENGTHSEQUENCEELEMENT_HPP
