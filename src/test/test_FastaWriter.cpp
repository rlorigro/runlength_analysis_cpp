#include "FastaWriter.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <assert.h>

using std::cout;
using std::experimental::filesystem::path;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_output_path = "/output/";
    path filename = "test_FastaWriter.fasta";
    path absolute_output_path = project_directory / relative_output_path / filename;

    cout << "WRITING: " << absolute_output_path << "\n";

    FastaWriter fasta_writer(absolute_output_path);
    sequence_element element;

    element = {.name="test1", .sequence="ACCAAACCCC"};
    fasta_writer.write(element);

    element = {.name="test2", .sequence="GGGGTTTGGT"};
    fasta_writer.write(element);
}
