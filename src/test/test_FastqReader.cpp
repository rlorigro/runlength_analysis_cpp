
#include "FastqReader.hpp"
#include <iostream>

using std::cout;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/test_sequences.fastq";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "Reading file: " << absolute_data_path << '\n';
    FastqReader reader(absolute_data_path);
    reader.filter_by_quality(7);

    return 0;
}
