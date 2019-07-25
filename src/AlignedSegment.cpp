
#include "AlignedSegment.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <string>

using std::string;

const array <string,10> Cigar::cigar_key = {"M","I","D","N","S","H","P","=","X"};

string Cigar::to_string() {
    return "(" + Cigar::cigar_key[this->code] + "," + std::to_string(this->length) + ")";
}

string AlignedSegment::to_string(){
    string reversal_string = reversal ? "R" : "F";

    string s = this->ref_name + " "
               + std::to_string(this->ref_start_index) + " "
               + reversal_string + " "
               + std::to_string(this->read_length) + " "
               + this->read_name + " "
               + std::to_string(this->n_cigar) + "\n";

    return s;
}