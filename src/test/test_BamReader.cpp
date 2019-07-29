#include "AlignedSegment.hpp"
#include "BamReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cassert>

using std::cout;
using std::experimental::filesystem::path;


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_data_path = "/data/test/test_alignable_sequences_non_RLE_VS_test_alignable_reference_non_RLE.sorted.bam";
    path absolute_data_path = project_directory / relative_data_path;

    cout << "TESTING " << absolute_data_path << "\n";
    cout << "BAM ITERATION TEST: \n";

    BamReader parser = BamReader(absolute_data_path);

    string region = "synthetic_ref_0";
    parser.initialize_region(region, 0, 1337);

    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    cout << cigar.to_string() << "\n";

    string ref_sequence = "ACCTTGCGATGCTAGCATAGCATGGCTCATGAATGCGATCCGATTGCAGTCCGAATGCATTGCTACGCAGTGCATAGGCTAGCTCAGACTGCTAGCTAGGCTACTAGCATGCCTAGTCAGTGACTAGCCTAGCGTTAGTCGATCATATCAGCGTACTCATCGATGCAGCACATGCATGCTATGTCTAGTACTACGGTACGATTATTCGATTCGCCGATGATCTAGCATGCGTACTGCTAGATGCTATGCATGCGTATGATATCTGATGCATGTCAGTTATGCATATTCGATATGTACTAGTTGCAGTCATGTGCATTATGCAGCTATTATTACGCTGAGTGCATAGCATGTCTGTCGCTAGCTAGATCGTAGCATGATCAGCATCATTCATGCATTCGTAGCATGCTAATGTCATTATAGCTATCGGCATTATCAGTAGCATCTAGCATAGATCTCAGTACGTAGTATCTATCAGTCAGTAGCTAGTCGATACGTATCGTTGCATATGCTACGATGATGCTATGCATGGCAATGCATTGCACTAGCTACGTATACTATGCATGTCAGTAGATGCTATACTGATCGAATCTTGATCTAGTAGCCTAGTAGCACTGACTGCATGTCAGTACGTAGCTACGTTCGTATACGATCATCTAGATCATGCTAGCATGCATGCATATATGTGACTGATGCTGATGCCGGATTATCGCGTATCGATCGATCATCATGATCATGATGATGCTGCATCAGAATACGCTGACTGACATCGACGACTGCATCGCGACTGCATCGGCAGCTAGCATGCGCGATGCATGCATCGTACTGCATGCAGTCGATCGATGCACGATGCATGCATGCATCGAGTATAGCCGGATTAGCTACTGAGCGATTTATCTCTGAGGAGATCTCGATCGTAGCATGCTGCGCATCTGCTAATGTCGGATGCTAGCGCTAGCTGCTTAGCTCATATTACGTATCTGATCTGATTCGATGCATGCATTATCGATTCGTATTAGCATCGTACGTAGCTATGCATTCGTAGCTAGCATCGTAGCTGAGCGATGCTATGCGCTAGCTTAGCGATGCTGCCGATCGTAGCGTATCAGAGTCGATCGTAGCTAGCTACGCGTACTAGCTAGCTACGTTAGCGCTAGATGATCTAGGCGCTATTATCGAGAGTCTCTAGGCTACTGATATCTGAGCAGGAAGAGTCGATCGTATGCTGCTGCTAGTCGTACGTATCGTATCGATGCATGTCATGCATAGTATGCGATCGCATGCTACTGTGCTGATGCTAGCTAGCTAGTCGATGTCGTAGCGGCATGTAGCGTACGCGG";

    string cigars;
    string ref_alignment;
    string read_alignment;

    while (parser.next_alignment(aligned_segment)){
        cigars = "";
        ref_alignment = "";
        read_alignment = "";

        int i = 0;

        cout << aligned_segment.to_string();

        aligned_segment.initialize_cigar_iterator();
        while (aligned_segment.next_coordinate_pair(coordinate, cigar)){
//            cout << cigar.to_string() << " " << coordinate.ref_index << " " << coordinate.read_index << "\n";
            cigars += cigar.get_cigar_code_as_string();

            if (cigar.is_ref_move()){
                ref_alignment += ref_sequence[coordinate.ref_index];
            }
            else{
                ref_alignment += "-";
            }

            if (cigar.is_read_move()){
                read_alignment += aligned_segment.get_read_base(coordinate.read_index);
            }
            else{
                read_alignment += "*";
            }

            i++;
        }

        cout << cigars << "\n";
        cout << ref_alignment << "\n";
        cout << read_alignment << "\n";
        cout << "\n";
    }
}
