#include "Identity.hpp"
#include "BedReader.hpp"

#include "boost/program_options.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"
#include <iostream>
#include <experimental/filesystem>
#include <utility>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <thread>
#include <exception>


using std::cout;
using std::cerr;
using std::flush;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::program_options::bool_switch;
using std::experimental::filesystem::path;
using boost::icl::interval_map;
using boost::icl::interval;
using boost::icl::partial_enricher;
using boost::icl::total_enricher;
using std::make_pair;
using std::vector;
using std::atomic;
using std::mutex;
using std::atomic_fetch_add;
using std::thread;
using std::move;
using std::ref;
using std::exception;
using std::to_string;



void parse_cigars(path bam_path,
                  unordered_map<string,FastaIndex>& ref_index,
                  regional_interval_map& bed_regions){

    // Initialize FastaReader and relevant containers
    SequenceElement sequence;

    // Initialize BAM reader and relevant containers
    BamReader bam_reader = BamReader(bam_path);
    AlignedSegment aligned_segment;
    Coordinate coordinate;
    Cigar cigar;

    string ref_name;

    bool filter_secondary = true;
    bool filter_supplementary = false;
    uint16_t map_quality_cutoff = 5;

    // Temporaries
    string true_base;
    string observed_base;

    uint8_t ambiguous_match_code = Cigar::cigar_code_key.at("M");
    uint8_t match_code = Cigar::cigar_code_key.at("=");
    uint8_t mismatch_code = Cigar::cigar_code_key.at("X");
    uint8_t insert_code = Cigar::cigar_code_key.at("I");
    uint8_t delete_code = Cigar::cigar_code_key.at("D");

    // Only allow matches
    unordered_set<uint8_t> valid_cigar_codes = {ambiguous_match_code,
                                                match_code,
                                                mismatch_code,
                                                insert_code,
                                                delete_code};
    string region_name;
    FastaIndex ref_fasta_index;

    for (auto& item: ref_index) {
        region_name = item.first;

        ref_fasta_index = item.second;

        // BAM coords are 1 based
        bam_reader.initialize_region(region_name, 1, ref_fasta_index.length);

        int i = 0;

        while (bam_reader.next_alignment(aligned_segment, map_quality_cutoff, filter_secondary, filter_supplementary)) {
            // Iterate cigars that match the criteria (must be '=')
            while (aligned_segment.next_coordinate(coordinate, cigar, valid_cigar_codes)) {

                if (bed_regions.count(region_name) == 0) {
                    continue;
                }

                auto result = bed_regions.at(region_name).find(coordinate.ref_index);

                if (coordinate.ref_index > 19000 and coordinate.ref_index < 19010){
                    cout << aligned_segment.read_name << '\t' << coordinate.read_true_index << '\n';
                }

                if (result != bed_regions.at(region_name).end()) {
                    cout << region_name << '\t' << result->first << '\t' << aligned_segment.read_name << '\t' << coordinate.read_true_index << '\n';
                }

                coordinate = {};
                cigar = {};
            }

            i++;
        }

        cerr << "\33[2K\rParsed: " << region_name << " " << to_string(ref_fasta_index.length) << flush;
    }
    cerr << '\n';
}


void extract_coordinates_from_alignment(path reference_fasta_path, path bam_path, path bed_path){
    // Initialize readers
    FastaReader ref_fasta_reader = FastaReader(reference_fasta_path);
    ref_fasta_reader.index();

    // Get reference contig lengths
    unordered_map<string,FastaIndex> ref_index;
    ref_index = ref_fasta_reader.get_index();

    // Parse the BED file as a map of interval maps (one interval map for each chromosome/contig/region)
    BedReader bed_reader(bed_path);
    regional_interval_map regions;
    bed_reader.read_regions(regions);


    for(auto [name, imap]: regions){
        for(auto interval: imap){
            cout << name << '\t' << interval.first << '\t' << interval.second << '\n';
        }
    }

    // Iterate all the alignments and locate read coordinates
    parse_cigars(bam_path, ref_index, regions);
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path bam_path;
    path bed_path;

    options_description options("Arguments");

    options.add_options()
            ("coordinates",
             value<path>(&bed_path),
             "File path of BED file containing coordinates to be found in the reads")

            ("bam",
             value<path>(&bam_path),
             "File path of BAM alignment file containing aligned sequences")

            ("ref",
             value<path>(&ref_fasta_path),
             "File path of FASTA file containing reference sequences");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    extract_coordinates_from_alignment(ref_fasta_path,
            bam_path,
            bed_path);

    return 0;
}

