#include "FastaWriter.hpp"
#include "htslib/hts.h"
//#include "hts.h"
#include "htslib/faidx.h"
//#include "faidx.h"
#include <string>
#include <iostream>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::cerr;
using std::experimental::filesystem::path;


int main(){
    // Find absolute path to test file within repo
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_input_path = "/data/test/";
    path filename = "test_sequences.fasta";
    path absolute_input_path = project_directory / relative_input_path / filename;

    // Build faidx
    int fai_exit_code = fai_build(absolute_input_path.c_str());
    if(fai_exit_code != 0) {
        cerr << "Error running faidx_build on " << absolute_input_path << "\n";
        exit(fai_exit_code);
    }

    // Load faidx
//    faidx_t m_fai = NULL;
//    m_fai = fai_load(absolute_input_path.c_str());

    // Fetch some stuff
//    auto& sequence = fai_fetch(m_fai, read_id.c_str(), &length);
//
//    cout << sequence;
}
