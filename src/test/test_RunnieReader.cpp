#include "DiscreteWeibull.hpp"
#include "RunnieReader.hpp"
#include "Matrix.hpp"
#include <iostream>

using std::cout;
using std::tie;


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/runnie/v2.1/";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "RUNNIE READER TEST: \n";

    RunnieReader reader = RunnieReader(absolute_data_path);
    reader.index();

    vector<RunnieSequenceElement> sequences;
    reader.fetch_all_sequences(sequences);

    runlength_matrix matrix(boost::extents[2][2][6][5]);

    for (auto& sequence: sequences) {
        cout << sequence.sequence << "\n";

        for (auto &scale: sequence.scales) {
            cout << scale << ",";
        }
        cout << "\n";

        for (auto &shape: sequence.shapes) {
            cout << shape << ",";
        }
        cout << "\n";

        vector<double> y(20);

        for (size_t i=0; i<sequence.shapes.size(); i++) {
            cout << "\n";
            evaluate_discrete_weibull(y, sequence.scales[i], sequence.shapes[i]);
            print_distribution(y);
        }

        bool reversal = true;
        uint8_t base_index = 0;
        uint16_t length = 5;
        update_runlength_matrix_with_weibull_probabilities(matrix, reversal, base_index, length, sequence.scales[0], sequence.shapes[0]);
        cout << matrix_to_string(matrix) << "\n";
    }

    vector<string> names = {"0e51e72e-4704-4629-a164-a2f2c4c77c07"};

    RunnieSequenceElement sequence;
    for (auto& name: names){
        cout << "TESTING " + name + ":\n";
        reader.fetch_sequence_bases(sequence, name);

        cout << sequence.sequence << "\n";
    }

    return 0;
}

