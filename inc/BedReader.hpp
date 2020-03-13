#ifndef RUNLENGTH_ANALYSIS_BEDREADER_H
#define RUNLENGTH_ANALYSIS_BEDREADER_H

#include "Region.hpp"
#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"

using boost::icl::interval_map;
using boost::icl::interval;
using boost::icl::total_enricher;
using std::experimental::filesystem::path;
using std::ifstream;
using std::string;
using std::vector;
using std::set;
using std::map;


using regional_interval_map = map <string, interval_map<uint64_t,bool,total_enricher> >;


class BedReader {
public:
    /// Attributes ///
    path bed_path;

    /// Methods ///
    BedReader(path bed_path);
    void read_regions(vector<Region>& regions);
    void read_regions(regional_interval_map& regions);
    bool next_line(Region& region);
    void subset_by_regions_name(vector<Region>& regions, set<string> names);

private:

    /// Attributes ///
    ifstream bed_file;
    static const size_t NAME = 0;
    static const size_t START = 1;
    static const size_t STOP = 2;
};


#endif //RUNLENGTH_ANALYSIS_BEDREADER_H
