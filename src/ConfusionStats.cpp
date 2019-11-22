#include "ConfusionStats.hpp"
#include "MarginPolishReader.hpp"
#include "ShastaReader.hpp"
#include "AlignedSegment.hpp"
#include "RunnieReader.hpp"
#include "FastaReader.hpp"
#include "FastaWriter.hpp"
#include "Runlength.hpp"
#include "BedReader.hpp"
#include "BamReader.hpp"
#include "Align.hpp"
#include "Base.hpp"
#include <vector>
#include <thread>
#include <string>
#include <iostream>
#include <mutex>
#include <exception>
#include <atomic>
#include <array>

using std::vector;
using std::string;
using std::thread;
using std::cout;
using std::cerr;
using std::flush;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::ref;
using std::move;
using std::exception;
using std::atomic;
using std::atomic_fetch_add;
using std::array;


void ConfusionStats::write_to_file(path output_file_path){
    cerr << "Writing file: " << output_file_path.string() << '\n';
    ofstream file(output_file_path);

    file << ">base_matches\n";
    for (auto& parent: this->base_match_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            file << index_to_base(parent_key) << ',' << child_key << ',' << child.second << '\n';
        }
    }

    file << ">base_mismatches\n";
    for (auto& parent: this->base_mismatch_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            file << index_to_base(parent_key) << ',' << child_key << ',' << child.second << '\n';
        }
    }

    file << ">length_matches\n";
    for (auto& parent: this->length_match_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            file << parent_key << ',' << child_key << ',' << child.second << '\n';
        }
    }

    file << ">length_mismatches\n";
    for (auto& parent: this->length_mismatch_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            file << parent_key << ',' << child_key << ',' << child.second << '\n';
        }
    }
}


uint16_t ConfusionStats::find_max_coverage(){
    uint16_t max_coverage = 0;
    uint16_t coverage = 0;

    for (auto& parent: this->base_mismatch_coverage) {
        coverage = parent.second.rbegin()->first;
        if (coverage > max_coverage){
            max_coverage = coverage;
        }
    }

    for (auto& parent: this->base_match_coverage) {
        coverage = parent.second.rbegin()->first;
        if (coverage > max_coverage){
            max_coverage = coverage;
        }
    }

    for (auto& parent: this->length_mismatch_coverage) {
        coverage = parent.second.rbegin()->first;
        if (coverage > max_coverage){
            max_coverage = coverage;
        }
    }

    for (auto& parent: this->length_match_coverage) {
        coverage = parent.second.rbegin()->first;
        if (coverage > max_coverage){
            max_coverage = coverage;
        }
    }

    return max_coverage;
}

void ConfusionStats::write_summary_to_file(path output_file_path){
    cerr << "Writing file: " << output_file_path.string() << '\n';
    ofstream file(output_file_path);

    const uint64_t max_coverage = this->find_max_coverage() + 1;
    vector<uint64_t > base_match_vector(max_coverage, 0);
    for (auto& parent: this->base_match_coverage){
        for (auto& child: parent.second){
            auto child_key = child.first;
            base_match_vector[child_key] += child.second;
        }
    }

    vector<uint64_t > base_mismatch_vector(max_coverage, 0);
    for (auto& parent: this->base_mismatch_coverage){
        for (auto& child: parent.second){
            auto child_key = child.first;
            base_mismatch_vector[child_key] += child.second;
        }
    }

    vector<uint64_t > length_match_vector(max_coverage, 0);
    for (auto& parent: this->length_match_coverage){
        for (auto& child: parent.second){
            auto child_key = child.first;
            length_match_vector[child_key] += child.second;
        }
    }

    vector<uint64_t > length_mismatch_vector(max_coverage, 0);
    for (auto& parent: this->length_mismatch_coverage){
        for (auto& child: parent.second){
            auto child_key = child.first;
            length_mismatch_vector[child_key] += child.second;
        }
    }


    file << ">base_matches\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        file << i << ',' << base_match_vector[i] << '\n';
    }
    file << ">base_mismatches\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        file << i << ',' << base_mismatch_vector[i] << '\n';
    }
    file << ">length_matches\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        file << i << ',' << length_match_vector[i] << '\n';
    }
    file << ">length_mismatches\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        file << i << ',' << length_mismatch_vector[i] << '\n';
    }

    float percent_error;

    file << ">base_error\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        percent_error = float(base_mismatch_vector[i]) / float(base_match_vector[i] + base_mismatch_vector[i]);
        file << i << ',' << percent_error << '\n';
    }
    file << ">length_error\n";
    for (uint64_t i=0; i<max_coverage; i++) {
        percent_error = float(length_mismatch_vector[i]) / float(length_match_vector[i] + length_mismatch_vector[i]);
        file << i << ',' << percent_error << '\n';
    }
}


void operator+=(ConfusionStats& a, ConfusionStats& b){

    for (auto& parent: b.base_match_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            a.base_match_coverage[parent_key][child_key] += child.second;
        }
    }

    for (auto& parent: b.base_mismatch_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            a.base_mismatch_coverage[parent_key][child_key] += child.second;
        }
    }

    for (auto& parent: b.length_match_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            a.length_match_coverage[parent_key][child_key] += child.second;
        }
    }

    for (auto& parent: b.length_mismatch_coverage){
        auto parent_key = parent.first;
        for (auto& child: parent.second){
            auto child_key = child.first;
            a.length_mismatch_coverage[parent_key][child_key] += child.second;
        }
    }
}


void ConfusionStats::update(
        char true_base,
        char consensus_base,
        uint16_t true_length,
        uint16_t consensus_length,
        uint16_t n_coverage){

    if (consensus_base == true_base){
        this->base_match_coverage[base_to_index(true_base)][n_coverage]++;
    }
    else {
        this->base_mismatch_coverage[base_to_index(true_base)][n_coverage]++;
    }

    if (consensus_length == true_length){
        this->length_match_coverage[true_length][n_coverage]++;
    }
    else {
        this->length_mismatch_coverage[true_length][n_coverage]++;
    }

}


template<typename T> void parse_aligned_coverage(path bam_path,
        path parent_directory,
        unordered_map <string,path>& read_paths,
        unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
        vector <Region>& regions,
        ConfusionStats& confusion_stats,
        atomic <uint64_t>& job_index){
    ///
    ///
    ///

    // Initialize Reader and relevant containers
    bool store_length_consensus = true;
    bool store_coverage_data = false;
    T reader = T(parent_directory, store_length_consensus, store_coverage_data);
    reader.set_index(read_paths);
    CoverageSegment segment;

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Region region;
    Cigar cigar;

    string ref_name;

    bool filter_secondary = true;
    uint16_t map_quality_cutoff = 5;

    bool in_left_bound;
    bool in_right_bound;
    char true_base;
    char consensus_base;
    uint16_t true_length = -1;
    uint16_t consensus_length = -1;
    uint16_t n_coverage;

    // Only allow matches and mismatches
    unordered_set<uint8_t> valid_cigar_codes = {Cigar::cigar_code_key.at("="),Cigar::cigar_code_key.at("X")};

    while (job_index < regions.size()) {
        uint64_t thread_job_index = job_index.fetch_add(1);
        region = regions.at(thread_job_index);

        // BAM coords are 1 based
        bam_reader.initialize_region(region.name, region.start+1, region.stop+1);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary)) {
            reader.fetch_read(segment, aligned_segment.read_name);

            // Iterate cigars that match the criteria (must be '=' or 'X')
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {
                in_left_bound = (int64_t(region.start) <= coordinate.ref_index - 1);
                in_right_bound = (coordinate.ref_index - 1 < int64_t(region.stop));

                // Subset alignment to portions of the read that are within the window/region
                if (in_left_bound and in_right_bound) {
                    true_base = ref_runlength_sequences.at(aligned_segment.ref_name).sequence[coordinate.ref_index];
                    consensus_base = segment.sequence[coordinate.read_true_index];

                    if (not is_valid_base(consensus_base) or not is_valid_base(true_base)){
                        continue;
                    }

                    if (aligned_segment.reversal){
                        consensus_base = complement_base(consensus_base);
                    }

                    true_length = ref_runlength_sequences.at(aligned_segment.ref_name).lengths[coordinate.ref_index];
                    consensus_length = segment.lengths[coordinate.read_true_index];

                    n_coverage = segment.n_coverage[coordinate.read_true_index];

                    confusion_stats.update(true_base, consensus_base, true_length, consensus_length, n_coverage);
                }
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region.to_string() << flush;
    }
}


template <typename T> ConfusionStats get_confusion_stats(path bam_path,
                                       path input_directory,
                                       unordered_map <string,path>& read_paths,
                                       unordered_map <string,RunlengthSequenceElement>& ref_runlength_sequences,
                                       vector <Region>& regions,
                                       uint16_t max_threads){
    ///
    ///
    ///

    vector<ConfusionStats> confusion_stats_per_thread(max_threads);

    vector<thread> threads;
    atomic<uint64_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<max_threads; i++){
        try {
            // Call thread safe function to read and write to file
            threads.emplace_back(thread(parse_aligned_coverage<T>,
                                        ref(bam_path),
                                        ref(input_directory),
                                        ref(read_paths),
                                        ref(ref_runlength_sequences),
                                        ref(regions),
                                        ref(confusion_stats_per_thread[i]),
                                        ref(job_index)));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }
    cerr << "\n" << flush;

    cerr << "Summing matrices from " << max_threads << " threads...\n";

    ConfusionStats confusion_sum;

    for (auto& item: confusion_stats_per_thread){
        confusion_sum += item;
    }

    confusion_stats_per_thread[0].write_summary_to_file("coverage_confusion_test/test_summary.csv");

    return confusion_sum;
}


void chunk_regions(path bed_path, vector<Region>& regions, unordered_map<string,RunlengthSequenceElement>& sequences, uint64_t chunk_size){
    if (bed_path.empty()){
        // Chunk alignment regions
        for (auto& [name, item]: sequences){
            chunk_sequence(regions, name, chunk_size, item.sequence.size());
        }
    }
    else{
        // Load the BED regions, and then find the union of the regions in BED and reference FASTA
        set<string> names;
        for (auto& item: sequences){
            names.insert(item.first);
        }

        BedReader bed_reader(bed_path);
        bed_reader.read_regions(regions);
        bed_reader.subset_by_regions_name(regions, names);
    }

}


template <typename T> void measure_confusion_stats_from_coverage_data(path input_directory,
                                                       path reference_fasta_path,
                                                       path output_directory,
                                                       uint16_t max_threads,
                                                       path bed_path){

    cerr << "Using " + to_string(max_threads) + " threads\n";

    // How big (bp) should the regions be for iterating the BAM? Regardless of size,
    // only one alignment worth of RAM is consumed per chunk. This value should be chosen as an appropriate
    // fraction of the genome size to prevent threads from being starved. Larger chunks also reduce overhead
    // associated with iterating reads that extend beyond the region (at the edges)
    uint64_t chunk_size = 1*1000*1000;

    T reader = T(input_directory);

    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);

    // Runlength encode the reference (store in memory)
    unordered_map<string,RunlengthSequenceElement> ref_runlength_sequences;
    path reference_fasta_path_rle;
    bool store_in_memory = true;
    reference_fasta_path_rle = runlength_encode_fasta_file(reference_fasta_path,
            ref_runlength_sequences,
            output_directory,
            store_in_memory,
            max_threads);

    reader.index();
    unordered_map<string,path> read_paths = reader.get_index();

    // Extract read names from index
    vector<string> read_names;
    for (auto& element: read_paths){
        read_names.push_back(element.first);
    }

    // Extract the sequences from MarginPolish TSVs (ignore coverage data for now)
    path reads_fasta_path_rle;
    reads_fasta_path_rle = write_all_consensus_sequences_to_fasta(reader,
            read_names,
            read_paths,
            input_directory,
            output_directory,
            max_threads);

    // Index all the files in the MP directory
    // Setup Alignment parameters
    bool sort = true;
    bool index = true;
    bool delete_intermediates = false;  //TODO: switch to true
    uint16_t k = 19;
    string minimap_preset = "asm20";    //TODO: make command line argument?
    bool explicit_mismatch = true;

    // Align coverage segments to the reference
    path bam_path;
    bam_path = align(reference_fasta_path_rle,
            reads_fasta_path_rle,
            output_directory,
            sort,
            index,
            delete_intermediates,
            k,
            minimap_preset,
            explicit_mismatch,
            max_threads);

    //TODO: cut off 50bp overlaps on MarginPolish TSVs?


    // If a BED file was provided, only iterate the regions of the BAM that may be found in the reference provided.
    // Otherwise, iterate the entire BAM.
    vector<Region> regions;
    chunk_regions(bed_path, regions, ref_runlength_sequences, chunk_size);

    cerr << "Iterating alignments...\n" << std::flush;

    // Launch threads for parsing alignments and generating matrices
    ConfusionStats stats = get_confusion_stats<T>(bam_path,
            input_directory,
            read_paths,
            ref_runlength_sequences,
            regions,
            max_threads);

    cerr << '\n';

    // Write output
    path output_file_path = output_directory / "confusion_stats.csv";
    output_file_path = absolute(output_file_path);
    stats.write_to_file(output_file_path);

    // Write output
    path summary_output_file_path = output_directory / "confusion_stats_summary.csv";
    summary_output_file_path = absolute(summary_output_file_path);
    stats.write_summary_to_file(summary_output_file_path);

}


void measure_confusion_stats_from_shasta(path input_directory,
        path reference_fasta_path,
        path output_directory,
        uint16_t max_threads,
        path bed_path) {

    measure_confusion_stats_from_coverage_data<ShastaReader>(input_directory,
            reference_fasta_path,
            output_directory,
            max_threads,
            bed_path);
}
