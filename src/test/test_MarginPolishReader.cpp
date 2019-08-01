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
    MarginPolishSegment segment;

    reader.index();
    map<string,path> read_paths = reader.get_index();
    string name;

    vector<MarginPolishSegment> segments;

    for (auto& element: read_paths){
        name = element.first;
        cout << element.first << " " << element.second << "\n";
        reader.fetch_read(segment, name);

        segment.print();
        segments.push_back(segment);
    }
}



