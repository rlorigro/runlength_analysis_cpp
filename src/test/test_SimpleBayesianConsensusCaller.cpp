
#include "SimpleBayesianConsensusCaller.hpp"
#include <iostream>
#include <vector>

using std::cout;
using std::vector;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path config_path = project_directory / "config/SimpleBayesianConsensusCaller-5.csv";

    SimpleBayesianConsensusCaller consensus_caller(config_path);
    vector <vector <float> > coverage;
    vector <float> consensus;

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {0, 1, 2},
            {3, 1, 2},
            {3, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {3, 1, 2},
            {3, 1, 3},
            {3, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {0, 1, 3},
            {0, 1, 3},
            {0, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {0, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {0, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    cout << "TEST\n";
    coverage = {
            {0, 1, 2},
            {0, 1, 2},
            {0, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
            {3, 1, 3},
    };

    consensus_caller(coverage, consensus);

    cout << "CONSENSUS: " << consensus[0] << " " << consensus[1] << '\n';

    return 0;
}

