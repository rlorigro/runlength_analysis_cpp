#include "Region.hpp"

Region::Region(string name, uint64_t start, uint64_t stop){
    this->name = name;
    this->start = start;
    this->stop = stop;
}


Region::Region() = default;


string Region::to_string(){
    string s = this->name + "_" + std::to_string(this->start) + "-" + std::to_string(this->stop);
    return s;
}
