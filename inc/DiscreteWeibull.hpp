
#ifndef RUNLENGTH_ANALYSIS_DISCRETEWEIBULL_H
#define RUNLENGTH_ANALYSIS_DISCRETEWEIBULL_H

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

using std::cout;
using std::string;
using std::vector;


double evaluate_weibull_cdf(size_t x, double scale, double shape);

void evaluate_discrete_weibull(vector<double>& distribution, double scale, double shape);

void print_distribution(vector<double>& distribution, uint16_t width=100, char character='|');


#endif //RUNLENGTH_ANALYSIS_DISCRETEWEIBULL_H
