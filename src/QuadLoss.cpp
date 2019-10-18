
#include "QuadLoss.hpp"

DiscreteWeibullLoss::DiscreteWeibullLoss(size_t distribution_size):
stats(distribution_size), distribution(distribution_size), distribution_size(distribution_size){}


DiscreteWeibullLoss::DiscreteWeibullLoss():
stats(0), distribution(0), distribution_size(0){}


void QuadLoss::reset(){
    return;
}

void QuadLoss::update(QuadCoordinate point){
    return;
}

double QuadLoss::calculate_loss(){
    return 0;
}

void DiscreteWeibullLoss::reset(){
    this->stats = MultiDistributionStats(this->distribution_size);
    this->distribution = vector<double>(this->distribution_size);
}

void DiscreteWeibullLoss::update(QuadCoordinate point){
    evaluate_discrete_weibull(this->distribution, point.x, point.y);
    this->stats.update(this->distribution);
}


double DiscreteWeibullLoss::calculate_loss(){
    if (this->stats.size == 0){
        return 0;
    }
    else {
//        return this->stats.calculate_loss();
        return this->stats.get_weighted_variance();
    }
}

