#include "MarginPolishReader.hpp"
#include <iostream>

using std::cout;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/marginpolish/";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";

    cout << "MARGINPOLISH READER TEST: \n";

    MarginPolishReader reader = MarginPolishReader(absolute_data_path);

    reader.index();
    map<string,path> read_paths = reader.get_index();

    for (auto& element: read_paths){
        cout << element.first << " " << element.second << "\n";
    }
}



