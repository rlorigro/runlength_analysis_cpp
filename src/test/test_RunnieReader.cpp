#include "RunnieReader.hpp"
#include <iostream>
#include <utility>

using std::cout;
using std::tie;


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/runnie/";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "RUNNIE READER TEST: \n";

    RunnieReader reader = RunnieReader(absolute_data_path);
    reader.index();
}

