#include <vector>
#include <cmath>

using std::vector;
using std::pow;
using std::exp;


double evaluate_weibull_cdf(size_t x, double scale, double shape){
    ///
    /// Evaluate the CDF of a given weibull distribution (defined by shape and scale) at a given value x
    ///

    return exp(-pow(x/scale, shape));
}


void evaluate_discrete_weibull(vector<double> distribution, double scale, double shape){
    ///
    /// Evaluate the Discrete Weibull Distribution for a given scale and shape, taking an ALLOCATED (!) vector which
    /// defines the range to evaluate over
    ///

    for (size_t x=0; x<distribution.size(); x++){
        distribution[x] = evaluate_weibull_cdf(x, scale, shape) - evaluate_weibull_cdf(x+1, scale, shape);
    }
}


