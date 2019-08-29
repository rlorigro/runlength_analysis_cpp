#include "FastaReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::experimental::filesystem::path;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/test_sequences.fasta";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";

    cout << "FASTA ITERATION TEST: \n";

    FastaReader fasta_reader(absolute_data_path);

    string line;
    SequenceElement element;

    size_t sequence_index = 0;
    while (!fasta_reader.end_of_file) {
        element = {};
        fasta_reader.next_element(element);

        if (sequence_index==0) {
            cout << "Testing first sequence: ";
            assert(element.name == "test1");
            assert(element.sequence == "ACCAAACCCC");
            cout << "PASS\n";
        }
        else if (sequence_index==1) {
            cout << "Testing second sequence: ";
            assert(element.name == "test2");
            assert(element.sequence == "GGGGTTTGGT");
            cout << "PASS\n";
        }
        else if (sequence_index==2) {
            cout << "Testing third sequence: ";
            assert(element.name == "test3");
            assert(element.sequence == "GGGGTTTGGTGGGGTTTGGTGGGGTTTGGT");
            cout << "PASS\n";
        }
        else if (sequence_index==3) {
            cout << "Testing fourth sequence: ";
            assert(element.name == "test4");
            assert(element.sequence == "ACCAAACCCC");
            cout << "PASS\n";
        }
        else if (sequence_index==4) {
            cout << "Testing NO empty sequences: \n";
            assert(sequence_index<2);
            cout << "PASS\n";
        }
        sequence_index++;
    }

    cout << "\nFASTA INDEXED FETCH TEST: \n";

    string sequence_name;

    cout << "Testing third sequence: ";
    element = {};
    sequence_name = "test3";
    fasta_reader.fetch_sequence(element, sequence_name);

    assert(element.name == "test3");
    assert(element.sequence == "GGGGTTTGGTGGGGTTTGGTGGGGTTTGGT");
    cout << "PASS\n";


    cout << "Testing fourth sequence: ";
    element = {};
    sequence_name = "test4";
    fasta_reader.fetch_sequence(element, sequence_name);

    assert(element.name == "test4");
    assert(element.sequence == "ACCAAACCCC");
    cout << "PASS\n";


    cout << "Testing first sequence: ";
    element = {};
    sequence_name = "test1";
    fasta_reader.fetch_sequence(element, sequence_name);

    assert(element.name == "test1");
    assert(element.sequence == "ACCAAACCCC");
    cout << "PASS\n";


    cout << "Testing second sequence: ";
    element = {};
    sequence_name = "test2";
    fasta_reader.fetch_sequence(element, sequence_name);

    assert(element.name == "test2");
    assert(element.sequence == "GGGGTTTGGT");
    cout << "PASS\n";

    return 0;
}

