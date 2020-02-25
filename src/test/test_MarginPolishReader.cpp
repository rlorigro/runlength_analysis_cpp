#include "MarginPolishReader.hpp"
#include <iostream>
#include <utility>

using std::cout;
using std::tie;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/marginpolish/";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "MARGINPOLISH READER TEST: \n";

//    MarginPolishReader reader = MarginPolishReader("/home/ryan/data/Nanopore/Human/marginpolish/easy-HG00733-RC");
    MarginPolishReader reader = MarginPolishReader(absolute_data_path);
    CoverageSegment segment;

    reader.index();
    unordered_map<string,path> read_paths = reader.get_index();
    string name;

    vector<CoverageSegment> segments;

    for (auto& element: read_paths){
        name = element.first;
        cout << "\n" << element.first << " " << element.second << "\n";
        reader.fetch_read(segment, name);

        segment.print();
        segments.push_back(segment);
    }

    for (auto& element: read_paths){
        name = element.first;
        cout << element.first << " " << element.second << "\n";
        reader.fetch_consensus_sequence(segment, name);

        segment.print();
        segments.push_back(segment);
    }


    path absolute_data_path_erroneous = project_directory / relative_data_path / "erroneous";

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "MARGINPOLISH READER ERROR INPUT TEST: \n";

//    MarginPolishReader reader = MarginPolishReader("/home/ryan/data/Nanopore/Human/marginpolish/easy-HG00733-RC");
    reader = MarginPolishReader(absolute_data_path_erroneous);
    segment = {};

    reader.index();
    read_paths = reader.get_index();
    segments.clear();

    for (auto& element: read_paths){
        name = element.first;
        cout << "\n" << element.first << " " << element.second << "\n";
        reader.fetch_read(segment, name);

        segment.print();
        segments.push_back(segment);
    }
    return 0;
}
