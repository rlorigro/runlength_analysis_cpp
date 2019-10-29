#include "RunlengthIndex.hpp"


ostream& operator<<(ostream& s, RunlengthIndex& index) {
    s << "sequence_length:\t" << index.sequence_length << '\n';
    s << "sequence_byte_index:\t" << index.sequence_byte_index << '\n';
    s << "name_length:\t" << index.name_length << '\n';
    s << "name:\t" << index.name;

    return s;
}

