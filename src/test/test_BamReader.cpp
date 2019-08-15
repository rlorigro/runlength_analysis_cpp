#include "AlignedSegment.hpp"
#include "BamReader.hpp"
#include "FastaReader.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <cassert>
#include <Runlength.hpp>

using std::cout;
using std::experimental::filesystem::path;


// Set of hand-curated cigar operations that should be observed in the test SAM
unordered_map <string, vector <pair <string, int64_t> > >  truth_set =
        {{"synthetic_read", {{"=",1337}}},
         {"synthetic_read_reverse", {{"=",1337}}},
         {"synthetic_read_delete_G_at_20", {{"=",19}, {"D",1}, {"=",1317}}},
         {"synthetic_read_delete_G_at_20_reverse", {{"=",19}, {"D",1}, {"=",1317}}},
         {"synthetic_read_insert_C_at_30", {{"=",29}, {"I",1}, {"=",1308}}},
         {"synthetic_read_insert_C_at_30_reverse", {{"=",29}, {"I",1}, {"=",1308}}},
         {"synthetic_read_sub_T_at_40", {{"=",40}, {"X",1}, {"=",1296}}},
         {"synthetic_read_sub_T_at_40_reverse", {{"=",40}, {"X",1}, {"=",1296}}},
         {"synthetic_read_N40", {{"S",41}, {"=",1296}}},
         {"synthetic_read_N40_reverse", {{"S",41}, {"=",1296}}},
         {"synthetic_read_500_clipped", {{"=",837}}},
         {"synthetic_read_500_clipped_reverse", {{"=",837}}},
         {"synthetic_read_300N_at_500", {{"=",499}, {"X",300}, {"=",538}}},
         {"synthetic_read_1500N_at_500", {{"S",1999}, {"=",538},{"=",499}, {"H",2038}}},            // Supplemental alignments combined
//         {"synthetic_read_1500N_at_500", {{"=",499}, {"H",2038}}},
         {"synthetic_read_1500N_at_500_reverse", {{"S",1999}, {"=",538},{"=",499}, {"H",2038}}},    // Supplemental alignments combined
//         {"synthetic_read_1500N_at_500_reverse", {{"=",499}, {"H",2038}}}
};


// Set of hand-curated cigar operations that should be observed in the test SAM
unordered_map <string, vector <pair <string, int64_t> > >  truth_set_sum =
        {{"synthetic_read", {{"=",1337}}},
         {"synthetic_read_reverse", {{"=",1337}}},
         {"synthetic_read_delete_G_at_20", {{"=",19+1317}, {"D",1}}},
         {"synthetic_read_delete_G_at_20_reverse", {{"=",19+1317}, {"D",1}}},
         {"synthetic_read_insert_C_at_30", {{"=",29+1308}, {"I",1}}},
         {"synthetic_read_insert_C_at_30_reverse", {{"=",29+1308}, {"I",1}}},
         {"synthetic_read_sub_T_at_40", {{"=",40+1296}, {"X",1}}},
         {"synthetic_read_sub_T_at_40_reverse", {{"=",40+1296}, {"X",1}}},
         {"synthetic_read_N40", {{"S",41}, {"=",1296}}},
         {"synthetic_read_N40_reverse", {{"S",41}, {"=",1296}}},
         {"synthetic_read_500_clipped", {{"=",837}}},
         {"synthetic_read_500_clipped_reverse", {{"=",837}}},
         {"synthetic_read_300N_at_500", {{"=",499+538}, {"X",300}}},
         {"synthetic_read_1500N_at_500", {{"=",538},{"=",499}}},            // Supplemental alignments combined
         {"synthetic_read_1500N_at_500_reverse", {{"=",538},{"=",499}}},            // Supplemental alignments combined
};


bool cigar_tuple_is_in_truth_set(unordered_map <string, vector <pair <string, int64_t> > >& truth, string read_name, string cigar_name, int64_t cigar_length){
    bool found = false;

    vector <pair <string, int64_t> > true_read_cigars = truth.at(read_name);

    for (const auto& element: true_read_cigars){
//        cout << element.first << element.second << "\n";
        if ((cigar_name == element.first) and (cigar_length == element.second)){
            found = true;
        }
    }

    return found;
}


bool test_cigars(unordered_map <string, vector <pair <string, int64_t> > >& truth, AlignedSegment& aligned_segment, string& cigars){
    SequenceElement cigar_sequence;
    cigar_sequence.name = aligned_segment.read_name;
    cigar_sequence.sequence = cigars;

    RunlengthSequenceElement cigar_sequence_RLE;
    runlength_encode(cigar_sequence_RLE, cigar_sequence);

//    cout << cigar_sequence_RLE.sequence << "\n";
//    for (auto& element: cigar_sequence_RLE.lengths){
//        cout << element << ",";
//    }
//    cout << "\n";

    for (size_t cigar_index = 0; cigar_index < cigar_sequence_RLE.sequence.size(); cigar_index++) {
        string cigar_name = string(1, cigar_sequence_RLE.sequence[cigar_index]);
        uint16_t cigar_length = cigar_sequence_RLE.lengths[cigar_index];
        bool cigar_found = cigar_tuple_is_in_truth_set(truth, cigar_sequence.name, cigar_name, cigar_length);

//        cout << cigar_index << "\t" << cigar_name << "\t" << cigar_length << "\t"
//             << cigar_found << "\n" << std::flush;

        assert(cigar_found);
    }

    cout << "PASS\n\n";
    return true;
}


int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test BAM path
    path relative_bam_path = "/data/test/test_alignable_sequences_non_RLE_VS_test_alignable_reference_non_RLE.sorted.bam";
    path absolute_bam_path = project_directory / relative_bam_path;

    // Get test FASTA reads path
    path relative_fasta_reads_path = "/data/test/test_alignable_sequences_non_RLE.fasta";
    path absolute_fasta_reads_path = project_directory / relative_fasta_reads_path;

    // Get test FASTA reference path
    path relative_fasta_ref_path = "/data/test/test_alignable_reference_non_RLE.fasta";
    path absolute_fasta_ref_path = project_directory / relative_fasta_ref_path;

    cout << "TESTING " << absolute_bam_path << "\n";
    cout << "BAM ITERATION TEST: \n";

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(absolute_bam_path);
    string ref_name = "synthetic_ref_0";
    bam_reader.initialize_region(ref_name, 0, 1337);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    // Initialize FASTA reader and relevant containers
    FastaReader reads_fasta_reader = FastaReader(absolute_fasta_reads_path);
    SequenceElement read_sequence;

    // Initialize FASTA reader and relevant containers
    FastaReader ref_fasta_reader = FastaReader(absolute_fasta_ref_path);
    SequenceElement ref_sequence;
    ref_fasta_reader.fetch_sequence(ref_sequence, ref_name);

    cout << cigar.to_string() << "\n";

    string cigars;
    string ref_alignment;
    string read_alignment;
    string read_alignment_inferred;

    while (bam_reader.next_alignment(aligned_segment)) {
        reads_fasta_reader.fetch_sequence(read_sequence, aligned_segment.read_name);

        cigars = "";
        ref_alignment = "";
        read_alignment = "";
        read_alignment_inferred = "";

        int i = 0;

        cout << aligned_segment.to_string() << "\n";

        while (aligned_segment.next_coordinate(coordinate, cigar)) {
            cigars += cigar.get_cigar_code_as_string();

            if (cigar.is_ref_move()) {
                ref_alignment += ref_sequence.sequence[coordinate.ref_index];
            } else {
                ref_alignment += "-";
            }
            if (cigar.is_read_move()) {
                read_alignment += aligned_segment.get_read_base(coordinate.read_index);
                read_alignment_inferred += read_sequence.sequence[coordinate.read_true_index];
            } else {
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

        test_cigars(truth_set, aligned_segment, cigars);
    }

    cout << "\n\nCIGAR SUBSET TEST:\n";

    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("M"),
//                                      Cigar::cigar_code_key.at("I"),
//                                      Cigar::cigar_code_key.at("D"),
//                                      Cigar::cigar_code_key.at("N"),
//                                      Cigar::cigar_code_key.at("S")
//                                      Cigar::cigar_code_key.at("H"),
//                                      Cigar::cigar_code_key.at("P"),
                                                Cigar::cigar_code_key.at("="),
//                                                Cigar::cigar_code_key.at("X")
    };

    bam_reader.initialize_region(ref_name, 0, 1337);

    while (bam_reader.next_alignment(aligned_segment)) {
        reads_fasta_reader.fetch_sequence(read_sequence, aligned_segment.read_name);

        cigars = "";
        ref_alignment = "";
        read_alignment = "";
        read_alignment_inferred = "";

        int i = 0;

        cout << aligned_segment.to_string() << "\n";

        while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
            if (i > 3000) break;
            cigars += cigar.get_cigar_code_as_string();

            if (cigar.is_ref_move()) {
                ref_alignment += ref_sequence.sequence[coordinate.ref_index];
            } else {
                ref_alignment += "-";
            }
            if (cigar.is_read_move()) {
                read_alignment += aligned_segment.get_read_base(coordinate.read_index);
                read_alignment_inferred += read_sequence.sequence[coordinate.read_true_index];
            } else {
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

        test_cigars(truth_set_sum, aligned_segment, cigars);
    }

    // Get test BAM path
    path relative_quality_bam_path = "/data/test/flag_and_quality_test.bam";
    path absolute_quality_bam_path = project_directory / relative_quality_bam_path;

    // Initialize BAM reader and relevant containers
    bam_reader = BamReader(absolute_quality_bam_path);
    ref_name = "synthetic_ref_0";
    bam_reader.initialize_region(ref_name, 0, 1337);

    uint16_t map_quality_threshold;
    bool filter_secondary;

    cout << "\nTESTING BAM FILTERS (NO FILTER):\n";

    map_quality_threshold = 0;
    filter_secondary = false;

    while (bam_reader.next_alignment(aligned_segment, map_quality_threshold, filter_secondary)) {
        cout << aligned_segment.to_string() << "\n";
    }

    cout << "\nTESTING BAM FILTERS (MQ>5 and NOT SECONDARY):\n";

    bam_reader.initialize_region(ref_name, 0, 1337);

    map_quality_threshold = 5;
    filter_secondary = true;

    while (bam_reader.next_alignment(aligned_segment, map_quality_threshold, filter_secondary)) {
        cout << aligned_segment.to_string() << "\n";
    }

    cout << "\nTESTING BAM FILTERS (MQ>5):\n";

    bam_reader.initialize_region(ref_name, 0, 1337);

    map_quality_threshold = 5;
    filter_secondary = false;

    while (bam_reader.next_alignment(aligned_segment, map_quality_threshold, filter_secondary)) {
        cout << aligned_segment.to_string() << "\n";
    }
}
