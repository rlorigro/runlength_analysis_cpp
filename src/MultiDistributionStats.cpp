
#include "MultiDistributionStats.hpp"
#include <cmath>

using std::pow;
using std::sqrt;


MultiDistributionStats::MultiDistributionStats():
    pointwise_stats(0){}


MultiDistributionStats::MultiDistributionStats(size_t distribution_size):
    pointwise_stats(distribution_size){}


void MultiDistributionStats::update(vector<double>& distribution){
    if (this->pointwise_stats.size() != distribution.size()){
        throw runtime_error("ERROR: user-provided distribution does not match MultiDistributionStats size");
    }

    for (size_t i=0; i<this->pointwise_stats.size(); i++){
        this->pointwise_stats[i].add(distribution[i]);
    }
}


double MultiDistributionStats::get_weighted_variance(){
    double variance = 0;
    double weight = 0;
    double sum = 0;

    for (auto& element: this->pointwise_stats){
        sum += element.get_mean();
    }

    for (auto& element: this->pointwise_stats){
        weight = element.get_mean()/sum;
        variance += element.get_variance()*weight;
    }

    return variance;
}


//double MultiDistributionStats::calculate_variance_from_pdf(){
//
//}


double MultiDistributionStats::get_meta_variance(){
    ///
    /// Variance = sum[ p(x)*(x-u)^2 ]
    ///
    size_t i;
    double sum = 0;
    double mean = 0;
    double variance = 0;
    double p = 0;

    for (auto& element: this->pointwise_stats){
        sum += element.get_mean();
    }

    i = 0;
    for (auto& element: this->pointwise_stats){
        p = element.get_mean() / sum;
        mean += p*i;

        i++;
    }

    i = 0;
    for (auto& element: this->pointwise_stats){
        p = element.get_mean() / sum;
        variance += p * pow(i-mean, 2);

        i++;
    }

    return variance;
}


double MultiDistributionStats::calculate_loss(){
    double noise = this->get_weighted_variance();
    double variance = this->get_meta_variance();

    double loss = sqrt(pow(noise,2) + pow(variance,2));

    return loss;
}


void MultiDistributionStats::get_pointwise_mean(vector<double>& means){
    for (size_t i=0; i<this->pointwise_stats.size(); i++){
        means[i] = this->pointwise_stats[i].get_mean();
    }
}


void MultiDistributionStats::get_pointwise_variance(vector<double>& variances){
    for (size_t i=0; i<this->pointwise_stats.size(); i++){
        variances[i] = this->pointwise_stats[i].get_variance();
    }
}
