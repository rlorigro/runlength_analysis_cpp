
#ifndef RUNLENGTH_ANALYSIS_RUNLENGTHINDEX_HPP
#define RUNLENGTH_ANALYSIS_RUNLENGTHINDEX_HPP

#include <experimental/filesystem>
#include <string>
#include <fstream>

using std::ostream;
using std::string;

class RunlengthIndex {
public:
    /// Attributes ///
    string name;
    uint64_t name_length;
    uint64_t sequence_byte_index;
    uint64_t sequence_length;

    /// Methods ///
};

ostream& operator<<(ostream& s, RunlengthIndex& index);

#endif //RUNLENGTH_ANALYSIS_RUNLENGTHINDEX_HPP
