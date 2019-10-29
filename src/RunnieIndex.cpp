#include "RunnieIndex.hpp"

using std::cout;

RunnieIndex::RunnieIndex(path file_path, uint64_t byte_index, uint64_t length){
    this->file_path = file_path;
    this->byte_index = byte_index;
    this->length = length;
}


ostream& RunnieIndex::operator<<(ostream& s){
    cout << this->file_path << " " << this->byte_index << " " << this->length;

    return s;
}

