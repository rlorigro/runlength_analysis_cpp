
#ifndef RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP
#define RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;


class SequenceElement{
public:
    string name;
    string sequence;

    vector<float> get_read_data();
    vector<float> generate_default_data_vector();
};


#endif //RUNLENGTH_ANALYSIS_SEQUENCEELEMENT_HPP
