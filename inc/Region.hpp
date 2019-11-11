#ifndef RUNLENGTH_ANALYSIS_REGION_HPP
#define RUNLENGTH_ANALYSIS_REGION_HPP

#include <string>

using std::string;


// For occasional convenience
class Region {
public:
    string name;
    uint64_t start;
    uint64_t stop;

    Region(string name, uint64_t start, uint64_t stop);
    Region();
    string to_string();
};


#endif //RUNLENGTH_ANALYSIS_REGION_HPP
