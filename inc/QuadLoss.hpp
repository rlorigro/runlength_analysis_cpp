
#ifndef RUNLENGTH_ANALYSIS_QUADLOSS_HPP
#define RUNLENGTH_ANALYSIS_QUADLOSS_HPP

#include "MultiDistributionStats.hpp"
#include "DiscreteWeibull.hpp"
#include "QuadTree.hpp"


class QuadLoss{
public:
    virtual void reset();
    virtual void update(QuadCoordinate point);
    virtual double calculate_loss();
};


class DiscreteWeibullLoss: public QuadLoss{
public:
    MultiDistributionStats stats;
    vector<double> distribution;
    size_t distribution_size;

    DiscreteWeibullLoss();
    DiscreteWeibullLoss(size_t distribution_size);
    void reset();
    void update(QuadCoordinate point);
    double calculate_loss();
};


#endif //RUNLENGTH_ANALYSIS_QUADLOSS_HPP
