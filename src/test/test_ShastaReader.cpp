#include "CoverageSegment.hpp"
#include "ShastaReader.hpp"
#include <iostream>
#include <utility>

using std::cout;
using std::tie;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/shasta/";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "SHASTA READER TEST: \n";

    ShastaReader reader = ShastaReader(absolute_data_path);
    CoverageSegment segment;

    reader.index();
    unordered_map<string,path> read_paths = reader.get_index();
    string name;

    vector<CoverageSegment> segments;

    for (auto& element: read_paths){
        name = element.first;
        cout << "\n" << name << " " << element.first << " " << element.second << "\n";
        reader.fetch_read(segment, name);

        segment.print();
        segments.push_back(segment);
    }

    for (auto& element: read_paths){
        name = element.first;
        cout << "\n" << name << " " << element.first << " " << element.second << "\n";
        reader.fetch_consensus_sequence(segment, name);

        segment.print();
        segments.push_back(segment);
    }

    return 0;
}
