
#include "SimpleBayesianRunnieConsensusCaller.hpp"
#include <iostream>
#include <vector>

using std::cout;
using std::vector;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-6-runnie-raw-reads-5mb-chr11.csv";

    SimpleBayesianRunnieConsensusCaller consensus_caller(config_path);
    vector <vector <float> > coverage;
    vector <float> consensus;

    cout << "TEST\n";
    coverage = {
            {0, 1, 2.1, 5.2},
            {0, 1, 2.1, 5.2},
            {0, 1, 2.1, 5.2},
            {3, 1, 2.1, 5.2},
            {3, 1, 3.1, 5.2},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';


    return 0;
}
