
#ifndef RUNLENGTH_ANALYSIS_MULTIDISTRIBUTIONSTATS_HPP
#define RUNLENGTH_ANALYSIS_MULTIDISTRIBUTIONSTATS_HPP

#include "IterativeSummaryStats.hpp"
#include <iostream>
#include <vector>
#include <cmath>

using std::vector;
using std::cout;

class MultiDistributionStats {
public:
    /// Attributes ///
    vector <IterativeSummaryStats<double> > pointwise_stats;

    /// Methods ///
    MultiDistributionStats();
    MultiDistributionStats(size_t distribution_size);
    void update(vector<double>& distribution);
    double get_weighted_variance();
    double get_meta_variance();
    double calculate_loss();
    void get_pointwise_mean(vector<double>& means);
    void get_pointwise_variance(vector<double>& variances);
};


#endif //RUNLENGTH_ANALYSIS_MULTIDISTRIBUTIONSTATS_HPP
