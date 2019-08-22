#include "DiscreteWeibull.hpp"
#include <cmath>

using std::pow;
using std::exp;
using std::to_string;


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

    for (size_t x=0; x<distribution.size()-1; x++){
        distribution[x+1] = evaluate_weibull_cdf(x, scale, shape) - evaluate_weibull_cdf(x+1, scale, shape);
    }
}


void print_distribution(vector<double>& distribution, uint16_t width, char character){
    ///
    /// Make a text representation of a distribution
    ///

    uint16_t n_characters;
    string line;

    size_t i = 0;
    for (auto& y: distribution){
        n_characters = uint16_t(y*width);
        line = to_string(i+1) + ":\t" + to_string(y) + "\t" + string(n_characters, character);
        cout << line << "\n";
        i++;
    }
}


