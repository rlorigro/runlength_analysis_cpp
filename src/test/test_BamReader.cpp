#include "BamReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::experimental::filesystem::path;


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/test_alignable_sequences_non_RLE_VS_test_alignable_reference_non_RLE.sorted.bam";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "FASTA ITERATION TEST: \n";

    BamReader parser = BamReader(absolute_data_path);

//    parser.initialize_region("synthetic_ref_0", 0, 1337);
    parser.print_region("synthetic_ref_0", 0, 1337);
}
