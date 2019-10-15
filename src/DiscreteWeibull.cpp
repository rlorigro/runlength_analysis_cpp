#include "DiscreteWeibull.hpp"
#include <stdexcept>
#include <cmath>

using std::pow;
using std::exp;
using std::to_string;
using std::runtime_error;


double evaluate_weibull_cdf(size_t x, double scale, double shape){
    ///
    /// Evaluate the CDF of a given weibull distribution (defined by shape and scale) at a given value x
    ///

    return exp(-pow(x/scale, shape));
}


void evaluate_discrete_weibull(vector<double>& distribution, double scale, double shape){
    ///
    /// Evaluate the Discrete Weibull Distribution for a given scale and shape, taking an ALLOCATED (!) vector which
    /// defines the range to evaluate over
    ///

    if (distribution.empty()){
        throw runtime_error("ERROR: uninitialized vector passed to 'evaluate_discrete_weibull', vector must be initialized");
    }

    for (size_t x=0; x<distribution.size()-1; x++){
        distribution[x+1] = evaluate_weibull_cdf(x, scale, shape) - evaluate_weibull_cdf(x+1, scale, shape);
    }
}


