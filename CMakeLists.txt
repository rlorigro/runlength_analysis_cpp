######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.0)
project(run_length_analysis VERSION 0.0.0)

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching -fsanitize=address")
set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------------------------------ #
##############################################

set(HTSLIB_HEADERS
        external/htslib
        external/htslib/htslib
        external/htslib/cram
        external/htslib/os
        )


set(HTSLIB_SOURCES
        external/htslib/cram/cram_codecs.c
        external/htslib/cram/cram_decode.c
        external/htslib/cram/cram_encode.c
        external/htslib/cram/cram_external.c
        external/htslib/cram/cram_index.c
        external/htslib/cram/cram_io.c
        external/htslib/cram/cram_samtools.c
        external/htslib/cram/cram_stats.c
        external/htslib/cram/files.c
        external/htslib/cram/mFILE.c
        external/htslib/cram/open_trace_file.c
        external/htslib/cram/pooled_alloc.c
        external/htslib/cram/rANS_static.c
        external/htslib/cram/sam_header.c
        external/htslib/cram/string_alloc.c
        external/htslib/os/rand.c
        external/htslib/test/fieldarith.c
        external/htslib/test/hfile.c
        external/htslib/test/hts_endian.c
        external/htslib/test/sam.c
        external/htslib/test/test-bcf-sr.c
        external/htslib/test/test-regidx.c
        external/htslib/test/test-vcf-api.c
        external/htslib/test/test-vcf-sweep.c
        external/htslib/test/test_bgzf.c
        external/htslib/test/test_view.c
        external/htslib/test/thrash_threads1.c
        external/htslib/test/thrash_threads2.c
        external/htslib/test/thrash_threads3.c
        external/htslib/test/thrash_threads4.c
        external/htslib/test/thrash_threads5.c
        external/htslib/test/thrash_threads6.c
        external/htslib/bcf_sr_sort.c
        external/htslib/bgzf.c
        external/htslib/bgzip.c
        external/htslib/errmod.c
        external/htslib/faidx.c
        external/htslib/hfile.c
        external/htslib/hfile_gcs.c
        external/htslib/hfile_libcurl.c
        external/htslib/hfile_net.c
        #external/htslib/hfile_s3.c
        external/htslib/hts.c
        external/htslib/htsfile.c
        external/htslib/kfunc.c
        external/htslib/knetfile.c
        external/htslib/kstring.c
        external/htslib/md5.c
        external/htslib/multipart.c
        #external/htslib/plugin.c
        external/htslib/probaln.c
        external/htslib/realn.c
        external/htslib/regidx.c
        external/htslib/sam.c
        external/htslib/synced_bcf_reader.c
        external/htslib/tabix.c
        external/htslib/tbx.c
        external/htslib/textutils.c
        external/htslib/thread_pool.c
        external/htslib/vcf.c
        external/htslib/vcf_sweep.c
        external/htslib/vcfutils.c
        )

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'curl' library")
endif()


include_directories(${HTSLIB_HEADERS})

# HTSLIB configure
add_custom_command(
        OUTPUT "${PROJECT_SOURCE_DIR}/external/htslib/config.h"
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-s3 --disable-plugins --disable-bz2
        COMMAND make
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/external/htslib/
)
add_custom_target(HTSLIB_CONFIGURED DEPENDS "${PROJECT_SOURCE_DIR}/external/htslib/config.h")
add_library(htslib ${HTSLIB_SOURCES})
add_dependencies(htslib HTSLIB_CONFIGURED)
target_link_libraries(htslib pthread bz2 z curl)



#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test executables.
set(SOURCES
        src/Align.cpp
        src/AlignedSegment.cpp
        src/BamReader.cpp
        src/Base.cpp
        src/BedReader.cpp
        src/BinaryIO.cpp
        src/BinaryRunnieWriter.cpp
        src/BinaryRunnieReader.cpp
        src/CigarKmer.cpp
        src/CompressedRunnieWriter.cpp
        src/CompressedRunnieReader.cpp
        src/ConfusionStats.cpp
        src/CoverageReader.cpp
        src/CoverageElement.cpp
        src/CoverageSegment.cpp
        src/DiscreteWeibull.cpp
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/FastqReader.cpp
        src/Identity.cpp
        src/IterativeSummaryStats.cpp
        src/Kmer.cpp
        src/MarginPolishReader.cpp
        src/Matrix.cpp
        src/Miscellaneous.cpp
        src/MultiDistributionStats.cpp
        src/Pileup.cpp
        src/PileupKmer.cpp
        src/PileupGenerator.cpp
        src/QuadCompressor.cpp
        src/QuadLoss.cpp
        src/QuadTree.cpp
        src/Quadrant.cpp
        src/Runlength.cpp
        src/RunnieIndex.cpp
        src/RunlengthSequenceElement.cpp
        src/RunnieSequenceElement.cpp
        src/ReferenceRunlength.cpp
        src/Region.cpp
        src/RunnieReader.cpp
        src/RunlengthWriter.cpp
        src/RunlengthReader.cpp
        src/RunlengthIndex.cpp
        src/SimpleBayesianConsensusCaller.cpp
        src/SimpleBayesianRunnieConsensusCaller.cpp
        src/SequenceElement.cpp
        src/ShastaReader.cpp
        )


project(runlength_analysis)
add_library(runlength_analysis SHARED ${SOURCES})

# To make sure the library is named runlength_analysis.so,
# get rid of the "lib" prefix.
set_target_properties(runlength_analysis PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(runlength_analysis PROPERTIES DEFINE_SYMBOL "")


###############################################
## ------------------------------------------ #
## -------- LINKING EXTERNAL LIBRARY -------- #
## ------------------------------------------ #
###############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_htslib
        URL https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
        CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_htslib/configure # --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/src/project_htslib-build/
        BUILD_IN_SOURCE True
        BUILD_COMMAND make all
        INSTALL_COMMAND ""
)

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_htslib INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/src/project_htslib/libhts.so)

# Define library as dependent on the downloaded project
add_dependencies(htslib project_htslib)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(runlength_analysis htslib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${INSTALL_DIR}/src/project_htslib/
        ${INSTALL_DIR}/src/project_htslib/os
        ${INSTALL_DIR}/src/project_htslib/cram
        ${INSTALL_DIR}/src/project_htslib/htslib)

## Add runtime path for main library so it can access external library
#set_property(TARGET runlength_analysis PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
##
##include_directories(${INSTALL_DIR}/src/project_htslib/)

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# stupid "experimental" filesystem library for gcc <8.0
link_libraries(stdc++fs)

find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# -------- TESTS --------

set(FILENAME_PREFIX test_htslib_indexed)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_FastaReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_FastaWriter)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_Runlength)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_BamReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_MarginPolishReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_Runlength_parallel)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_RunnieReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_PileupGenerator)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_boost_interval_map)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_boost_interval_map_bug)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_CompressedRunnieWriter)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_ShastaReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_ReferenceRunlength)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_QuadTree)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_IterativeSummaryStats)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_MultiDistributionStats)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_RunlengthWriter)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_BinaryRunnieWriter)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_SimpleBayesianConsensusCaller)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_predict_minimap_runnie_consensus_v2)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_predict_minimap_runlength_consensus_v2)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_SimpleBayesianRunnieConsensusCaller)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_logsumexp)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_BedReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_kmer)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_kmer_iterator)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_FastqReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# -------- SCRIPTS --------

set(FILENAME_PREFIX fasta_to_RLE_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX align)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_distribution_from_marginpolish)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_distribution_from_shasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_distribution_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_distribution_from_runnie)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX count_marginpolish_coverage)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX compress_runnie)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_identity_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_identity_per_alignment_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_identity_from_bam)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX load_runnie_as_quad_tree)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX quad_compress_runnie)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_confusion_stats_from_shasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_read_length_stats)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_kmer_coverage_stats_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_kmer_confusion_stats_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_runlength_kmer_confusion_stats_from_bam)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_kmer_confusion_stats_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX measure_read_kmer_identity_from_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX label_coverage_data_from_shasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX predict_minimap_runlength_consensus_by_coverage)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX predict_consensus_from_shasta_coverage)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX filter_fastq_by_quality)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

#set(FILENAME_PREFIX extract_read_coordinates_from_alignment)
#add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
#set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
#target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        runlength_analysis
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

