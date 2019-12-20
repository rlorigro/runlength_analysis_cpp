#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <experimental/filesystem>
#include "Miscellaneous.hpp"

using std::string;
using std::to_string;
using std::vector;
using std::remove;
using std::replace;
using std::cout;
using std::cerr;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::absolute;
using std::experimental::filesystem::create_directories;


path minimap_align(path ref_sequence_path,
                   path read_sequence_path,
                   path output_dir,
                   string minimap_preset,
                   bool explicit_mismatch,
                   uint16_t max_threads,
                   uint16_t k){

    // Find filename prefixes to be combined to generate predictable output filename
    string ref_filename_prefix;
    string read_filename_prefix;

    // This works because etc_prefix is a string object, which means the value is copied
    ref_filename_prefix = ref_sequence_path.filename().replace_extension("").string();
    replace(ref_filename_prefix.begin(), ref_filename_prefix.end(), '.', '_');

    // This works because etc_prefix is a string object, which means the value is copied
    read_filename_prefix = read_sequence_path.filename().replace_extension("").string();
    replace(read_filename_prefix.begin(), read_filename_prefix.end(), '.', '_');

    path output_filename = read_filename_prefix + "_VS_" + ref_filename_prefix + ".sam";
    path output_path = output_dir / output_filename;

    cerr << "REDIRECTING TO: " << output_filename.string() << "\n";

    vector<string> arguments;
    if (not minimap_preset.empty()) {
        // Set up arguments in a readable, modular format
        arguments = {"minimap2",
                                    "-a",
                                    "-t", to_string(max_threads),
                                    "-x", minimap_preset,
                                    "-k", to_string(k),
                                    ref_sequence_path.string(),
                                    read_sequence_path.string(),
                                    ">", output_path.string()};
    }
    else if (k==0){
        arguments = {"minimap2",
                                    "-a",
                                    "-t", to_string(max_threads),
                                    "-x", minimap_preset,
                                    ref_sequence_path.string(),
                                    read_sequence_path.string(),
                                    ">", output_path.string()};
    }
    else{
        arguments = {"minimap2",
                                    "-a",
                                    "-t", to_string(max_threads),
                                    "-k", to_string(k),
                                    ref_sequence_path.string(),
                                    read_sequence_path.string(),
                                    ">", output_path.string()};
    }

    if (explicit_mismatch){
        arguments.insert(arguments.begin() + 1, "--eqx");
    }

    // Convert arguments to single string
    string argument_string = join(arguments, ' ');
    cerr << "\nRUNNING: " << argument_string << "\n";

    run_command(argument_string);

    return output_path.string();
}


path samtools_sort(path input_path, uint16_t max_threads){
    path output_path = input_path;
    output_path.replace_extension(".sorted.bam");

    // Set up arguments in a readable, modular format
    vector <string> arguments = {"samtools sort",
                                 "-O", "BAM",
                                 "-@", to_string(max_threads),
                                 "-o", output_path.string(),
                                 input_path.string()};


    // Convert arguments to single string
    string argument_string = join(arguments, ' ');
    cerr << "\nRUNNING: " << argument_string << "\n";

    run_command(argument_string);

    return output_path.string();
}


path samtools_index(path input_path, uint16_t max_threads){
    path output_path = input_path.string() + ".bai";

    // Set up arguments in a readable, modular format
    vector <string> arguments = {"samtools index",
                                 "-@", to_string(max_threads),
                                 input_path};

    // Convert arguments to single string
    string argument_string = join(arguments, ' ');
    cerr << "\nRUNNING: " << argument_string << "\n";

    run_command(argument_string);

    return output_path.string();
}


path align(path ref_sequence_path,
           path read_sequence_path,
           path output_dir,
           bool sort,
           bool index,
           bool delete_intermediates,
           uint16_t k,
           string minimap_preset,
           bool explicit_mismatch,
           uint16_t max_threads){

    // Convert output dir to absolute dir
    output_dir = absolute(output_dir);

    // Ensure output dir exists
    create_directories(output_dir);

    // Create placeholders for all possible file intermediates
    path sam_output_path;
    path sorted_bam_output_path;
    path bai_output_path;
    path output_path;

    // Perform the initial call to aligner
    sam_output_path = minimap_align(ref_sequence_path, read_sequence_path, output_dir, minimap_preset, explicit_mismatch, max_threads, k);
    output_path = sam_output_path;

    if (sort) {
        // Do samtools sort
        sorted_bam_output_path = samtools_sort(sam_output_path, max_threads);
        output_path = sorted_bam_output_path;

        if (index) {
            // Do samtools index
            bai_output_path = samtools_index(sorted_bam_output_path, max_threads);
        }
        if (delete_intermediates) {
            // Delete uncompressed SAM
            remove(sam_output_path);
        }
    }

    // return the path to whichever sam/bam file was created last
    return output_path.string();
}
