#include "BedReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::experimental::filesystem::path;
using std::runtime_error;


int  main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/bed/test.bed";
    path absolute_data_path = project_directory / relative_data_path;

    BedReader reader(absolute_data_path);

    vector<Region> regions;
    reader.read_regions(regions);

    if (regions.size() != 4){
        throw runtime_error("FAIL: incorrect number of regions");
    }

    vector<Region> truth_set = {{"chr1", 109740, 109773},
                                {"chr1", 109893, 110138},
                                {"chr2", 12342, 12526},
                                {"chr2", 12629, 12655}};

    for (size_t i=0; i<regions.size(); i++){
        cout << regions[i].to_string() << '\n';
        cout << truth_set[i].to_string() << '\n';

        bool a = regions[i].name == truth_set[i].name;
        bool b = regions[i].start == truth_set[i].start;
        bool c = regions[i].stop == truth_set[i].stop;

        if (not (a and b and c)){
            throw runtime_error("FAIL: incorrect value in BED region");
        }
    }

    reader.subset_by_regions_name(regions, {"chr1"});

    truth_set = {{"chr1", 109740, 109773},
                 {"chr1", 109893, 110138}};

    for (size_t i=0; i<regions.size(); i++){
        bool a = regions[i].name == truth_set[i].name;
        bool b = regions[i].start == truth_set[i].start;
        bool c = regions[i].stop == truth_set[i].stop;

        if (not (a and b and c)){
            throw runtime_error("FAIL: incorrect value in BED region");
        }
    }
    cerr << "PASS\n";

    return 0;
}
