#include "PileupGenerator.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Runlength.hpp"
#include "Align.hpp"
#include "Kmer.hpp"
#include "Base.hpp"

using std::cerr;
using std::ifstream;


void generate_test_pileup(Pileup& pileup){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test FASTA reads path
    path relative_fasta_reads_path = "/data/test/test_alignable_sequences_non_RLE.fasta";
    path absolute_fasta_reads_path = project_directory / relative_fasta_reads_path;

    // Get test FASTA reference path
    path relative_fasta_ref_path = "/data/test/test_alignable_reference_non_RLE.fasta";
    path absolute_fasta_ref_path = project_directory / relative_fasta_ref_path;

    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;
    uint16_t k = 19;
    string minimap_preset = "asm20";
    bool explicit_mismatch = true;
    path output_directory = project_directory / "/data/test/";
    uint16_t max_threads = 30;

    // Align reads to the reference
    path bam_path;
//    bam_path = align(
//            absolute_fasta_ref_path,
//            absolute_fasta_reads_path,
//            output_directory,
//            sort,
//            index,
//            delete_intermediates,
//            k,
//            minimap_preset,
//            explicit_mismatch,
//            max_threads);

    bam_path = align(absolute_fasta_ref_path,
           absolute_fasta_reads_path,
           output_directory,
           sort,
           index,
           delete_intermediates,
           k,
           minimap_preset,
           explicit_mismatch,
           max_threads);

    path absolute_bam_path = absolute(bam_path);

    PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, 20);

    Region region = Region("synthetic_ref_0", 0, 1337);

    FastaReader ref_reader = FastaReader(absolute_fasta_ref_path);
    FastaReader sequence_reader = FastaReader(absolute_fasta_reads_path);

    pileup_generator.fetch_region(region, sequence_reader, pileup);
}


int main(){
    Pileup pileup;

    generate_test_pileup(pileup);

    PileupGenerator::print(pileup);

    size_t depth_index = 1;
    size_t window_size = 7;
    uint8_t k = 10;

    PileupReadKmerIterator iterator(window_size, depth_index, k);

    deque<uint8_t> kmer;

    // Test a single read iterator in the pileup
    for (size_t i=0; i<40 - 1; i++) {
        unordered_map <uint64_t, unordered_set <uint64_t> > middle_kmers;

        // Walk along the pileup in a windowed manner
        iterator.step(pileup, middle_kmers);

        for (const auto& item: iterator.current_kmer) {
            cout << index_to_base(item);
        }
        cout << '\n';

        for (const auto& [middle_kmer_index, coverage_set]: middle_kmers) {
            index_to_kmer(kmer, middle_kmer_index, k);
            for (auto& b: kmer) {
                cout << index_to_base(b);
            }
            cout << ':' << coverage_set.size() << ',';
        }
        cout << '\n' << '\n';
    }


    // Test the overarching class which contains as many read iterators as the max observed coverage in the pileup
    PileupKmerIterator pileup_iterator(pileup, window_size, k);

    cout << "max_depth: " << pileup.max_observed_depth << '\n';

    for (size_t i=0; i<pileup.pileup.size() - 1; i++) {
        pileup_iterator.step(pileup);
    }

    cout << pileup_iterator.kmer_stats.to_string();

    return 0;
}

