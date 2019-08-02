#include "AlignedSegment.hpp"
#include "BamReader.hpp"
#include "FastaReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cassert>

using std::cout;
using std::experimental::filesystem::path;


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test BAM path
    path relative_bam_path = "/data/test/test_alignable_sequences_non_RLE_VS_test_alignable_reference_non_RLE.sorted.bam";
    path absolute_bam_path = project_directory / relative_bam_path;

    // Get test FASTA path
    path relative_fasta_path = "/data/test/test_alignable_sequences_non_RLE.fasta";
    path absolute_fasta_path = project_directory / relative_fasta_path;

    cout << "TESTING " << absolute_bam_path << "\n";
    cout << "BAM ITERATION TEST: \n";

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(absolute_bam_path);
    string region = "synthetic_ref_0";
    bam_reader.initialize_region(region, 0, 1337);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    // Initialize FASTA reader and relevant containers
    FastaReader fasta_reader = FastaReader(absolute_fasta_path);
    SequenceElement fasta_sequence;

    cout << cigar.to_string() << "\n";

    string ref_sequence = "ACCTTGCGATGCTAGCATAGCATGGCTCATGAATGCGATCCGATTGCAGTCCGAATGCATTGCTACGCAGTGCATAGGCTAGCTCAGACTGCTAGCTAGGCTACTAGCATGCCTAGTCAGTGACTAGCCTAGCGTTAGTCGATCATATCAGCGTACTCATCGATGCAGCACATGCATGCTATGTCTAGTACTACGGTACGATTATTCGATTCGCCGATGATCTAGCATGCGTACTGCTAGATGCTATGCATGCGTATGATATCTGATGCATGTCAGTTATGCATATTCGATATGTACTAGTTGCAGTCATGTGCATTATGCAGCTATTATTACGCTGAGTGCATAGCATGTCTGTCGCTAGCTAGATCGTAGCATGATCAGCATCATTCATGCATTCGTAGCATGCTAATGTCATTATAGCTATCGGCATTATCAGTAGCATCTAGCATAGATCTCAGTACGTAGTATCTATCAGTCAGTAGCTAGTCGATACGTATCGTTGCATATGCTACGATGATGCTATGCATGGCAATGCATTGCACTAGCTACGTATACTATGCATGTCAGTAGATGCTATACTGATCGAATCTTGATCTAGTAGCCTAGTAGCACTGACTGCATGTCAGTACGTAGCTACGTTCGTATACGATCATCTAGATCATGCTAGCATGCATGCATATATGTGACTGATGCTGATGCCGGATTATCGCGTATCGATCGATCATCATGATCATGATGATGCTGCATCAGAATACGCTGACTGACATCGACGACTGCATCGCGACTGCATCGGCAGCTAGCATGCGCGATGCATGCATCGTACTGCATGCAGTCGATCGATGCACGATGCATGCATGCATCGAGTATAGCCGGATTAGCTACTGAGCGATTTATCTCTGAGGAGATCTCGATCGTAGCATGCTGCGCATCTGCTAATGTCGGATGCTAGCGCTAGCTGCTTAGCTCATATTACGTATCTGATCTGATTCGATGCATGCATTATCGATTCGTATTAGCATCGTACGTAGCTATGCATTCGTAGCTAGCATCGTAGCTGAGCGATGCTATGCGCTAGCTTAGCGATGCTGCCGATCGTAGCGTATCAGAGTCGATCGTAGCTAGCTACGCGTACTAGCTAGCTACGTTAGCGCTAGATGATCTAGGCGCTATTATCGAGAGTCTCTAGGCTACTGATATCTGAGCAGGAAGAGTCGATCGTATGCTGCTGCTAGTCGTACGTATCGTATCGATGCATGTCATGCATAGTATGCGATCGCATGCTACTGTGCTGATGCTAGCTAGCTAGTCGATGTCGTAGCGGCATGTAGCGTACGCGG";

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;

    while (bam_reader.next_alignment(aligned_segment)){
        fasta_reader.fetch_sequence(fasta_sequence, aligned_segment.read_name);

        cigars = "";
        ref_alignment = "";
        read_alignment = "";
        read_alignment_inferred = "";

        int i = 0;

        cout << aligned_segment.to_string();

        aligned_segment.initialize_cigar_iterator();
        while (aligned_segment.next_coordinate(coordinate, cigar)){
            cigars += cigar.get_cigar_code_as_string();

            if (cigar.is_ref_move()){
                ref_alignment += ref_sequence[coordinate.ref_index];
            }
            else{
                ref_alignment += "-";
            }

            if (cigar.is_read_move()){
                read_alignment += aligned_segment.get_read_base(coordinate.read_index);
                read_alignment_inferred += fasta_sequence.sequence[coordinate.read_true_index];
            }
            else{
                read_alignment += "*";
                read_alignment_inferred += "*";
            }

            i++;
        }

        cout << cigars << "\n";
        cout << ref_alignment << "\n";
        cout << read_alignment << "\n";
        cout << read_alignment_inferred << "\n";
        cout << "\n";
    }
}
