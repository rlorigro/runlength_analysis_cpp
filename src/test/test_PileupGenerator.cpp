
#include "PileupGenerator.hpp"


int main(){
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

    PileupGenerator pileup_generator = PileupGenerator(absolute_bam_path, absolute_fasta_ref_path, absolute_fasta_reads_path, 10);

    Region region = Region("synthetic_ref_0", 0, 1337);
    pileup_generator.fetch_region(region);

//    region = Region("synthetic_ref_0", 0, 1337);
//    pileup_generator.fetch_region(region, 10);

}
