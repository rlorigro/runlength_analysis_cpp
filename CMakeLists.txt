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
add_definitions(-ggdb3 -O0 -Wall)

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
        externalTools/htslib
        externalTools/htslib/htslib
        externalTools/htslib/cram
        externalTools/htslib/os
        )


set(HTSLIB_SOURCES
        externalTools/htslib/cram/cram_codecs.c
        externalTools/htslib/cram/cram_decode.c
        externalTools/htslib/cram/cram_encode.c
        externalTools/htslib/cram/cram_external.c
        externalTools/htslib/cram/cram_index.c
        externalTools/htslib/cram/cram_io.c
        externalTools/htslib/cram/cram_samtools.c
        externalTools/htslib/cram/cram_stats.c
        externalTools/htslib/cram/files.c
        externalTools/htslib/cram/mFILE.c
        externalTools/htslib/cram/open_trace_file.c
        externalTools/htslib/cram/pooled_alloc.c
        externalTools/htslib/cram/rANS_static.c
        externalTools/htslib/cram/sam_header.c
        externalTools/htslib/cram/string_alloc.c
        externalTools/htslib/os/rand.c
        externalTools/htslib/test/fieldarith.c
        externalTools/htslib/test/hfile.c
        externalTools/htslib/test/hts_endian.c
        externalTools/htslib/test/sam.c
        externalTools/htslib/test/test-bcf-sr.c
        externalTools/htslib/test/test-regidx.c
        externalTools/htslib/test/test-vcf-api.c
        externalTools/htslib/test/test-vcf-sweep.c
        externalTools/htslib/test/test_bgzf.c
        externalTools/htslib/test/test_view.c
        externalTools/htslib/test/thrash_threads1.c
        externalTools/htslib/test/thrash_threads2.c
        externalTools/htslib/test/thrash_threads3.c
        externalTools/htslib/test/thrash_threads4.c
        externalTools/htslib/test/thrash_threads5.c
        externalTools/htslib/test/thrash_threads6.c
        externalTools/htslib/bcf_sr_sort.c
        externalTools/htslib/bgzf.c
        externalTools/htslib/bgzip.c
        externalTools/htslib/errmod.c
        externalTools/htslib/faidx.c
        externalTools/htslib/hfile.c
        externalTools/htslib/hfile_gcs.c
        externalTools/htslib/hfile_libcurl.c
        externalTools/htslib/hfile_net.c
        #externalTools/htslib/hfile_s3.c
        externalTools/htslib/hts.c
        externalTools/htslib/htsfile.c
        externalTools/htslib/kfunc.c
        externalTools/htslib/knetfile.c
        externalTools/htslib/kstring.c
        externalTools/htslib/md5.c
        externalTools/htslib/multipart.c
        #externalTools/htslib/plugin.c
        externalTools/htslib/probaln.c
        externalTools/htslib/realn.c
        externalTools/htslib/regidx.c
        externalTools/htslib/sam.c
        externalTools/htslib/synced_bcf_reader.c
        externalTools/htslib/tabix.c
        externalTools/htslib/tbx.c
        externalTools/htslib/textutils.c
        externalTools/htslib/thread_pool.c
        externalTools/htslib/vcf.c
        externalTools/htslib/vcf_sweep.c
        externalTools/htslib/vcfutils.c
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
        OUTPUT "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h"
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-s3 --disable-plugins --disable-bz2
        COMMAND make
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/externalTools/htslib/
)
add_custom_target(HTSLIB_CONFIGURED DEPENDS "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h")
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
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/Runlength.cpp
        src/Miscellaneous.cpp
        src/Align.cpp
        src/AlignedSegment.cpp)
#        src/CigarParser.cpp)   # !!!!!!!!!!!!!!!!! UNCOMMENT THIS TO SEE HOW HTSLIB IS BROKEN !!!!!!!!!!!!!!!

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
#
## Need to explicitly enable ExternalProject functionality
#include(ExternalProject)
#
## Download or update library as an external project
#ExternalProject_Add(project_htslib
#        URL https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
#        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/
#        CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/external/src/project_htslib/configure # --prefix=${CMAKE_CURRENT_SOURCE_DIR}/external/src/project_htslib-build/
#        BUILD_IN_SOURCE True
#        BUILD_COMMAND make all
#        INSTALL_COMMAND ""
#)
#
## Define INSTALL_DIR as the install directory for external library
#ExternalProject_Get_Property(project_htslib INSTALL_DIR)
#
## Create new library for external project (so it can be linked with main library)
#add_library(htslib STATIC IMPORTED)
#set_property(TARGET htslib
#        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/src/project_htslib/libhts.so)
#
## Define library as dependent on the downloaded project
#add_dependencies(htslib project_htslib)
#
## Define main library as dependent on the downloaded project (transitively)
#add_dependencies(runlength_analysis htslib)
#
## Ensure that main library has access to primary dependencies' and secondary dependencies' headers
#include_directories(
#        ${INSTALL_DIR}/src/project_htslib/
#        ${INSTALL_DIR}/src/project_htslib/os
#        ${INSTALL_DIR}/src/project_htslib/cram
#        ${INSTALL_DIR}/src/project_htslib/htslib)
#
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

set(FILENAME_PREFIX test_htslib)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_htslib_indexed)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_threading_simple)
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

set(FILENAME_PREFIX test_fai_fetch)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX test_fai_fetch_custom)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# !!!!!!!!!!!!!  UNCOMMENT THIS TO SEE HOW HTSLIB IS BROKEN !!!!!!!!!!!!!!!!!!!!!!!
set(FILENAME_PREFIX test_CigarParser)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)


# -------- SCRIPTS --------

set(FILENAME_PREFIX fasta_to_RLE_fasta)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX fasta_to_RLE_fasta_parallel)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

set(FILENAME_PREFIX align)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/src/project_htslib/")
target_link_libraries(${FILENAME_PREFIX} runlength_analysis htslib Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        runlength_analysis
        test_htslib
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

