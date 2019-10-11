
#ifndef RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP
#define RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP

#include <iostream>
#include <stdexcept>

using std::runtime_error;


template <class T> class IterativeSummaryStats {
private:
    /// Attributes ///
    T k;
    T n;
    T sum;
    T sum_of_squares;

public:
    /// Methods ///
    IterativeSummaryStats();
    void add(T x);
    void remove(T x);
    double get_mean();
    double get_variance();
};

template <class T> IterativeSummaryStats<T>::IterativeSummaryStats(){
    this->k = 0;
    this->n = 0;
    this->sum = 0;
    this->sum_of_squares = 0;
}

template <class T> void IterativeSummaryStats<T>::add(T x){
    if (this->n == 0){
        this->k = x;
    }

    this->sum += x - this->k;
    this->sum_of_squares += (x - this->k)*(x - this->k);
    this->n++;
}

template <class T> void IterativeSummaryStats<T>::remove(T x){
    if (this->n == 0){
        throw runtime_error("ERROR: cannot remove a value from empty IterativeSummaryStats object");
    }

    this->sum -= x - this->k;
    this->sum_of_squares -= (x - this->k)*(x - this->k);
    this->n--;
}

template <class T> double IterativeSummaryStats<T>::get_mean() {
    ///
    /// K + Ex / n
    ///
    double mean = this->k + double(this->sum)/(this->n);

    return mean;
}

template <class T> double IterativeSummaryStats<T>::get_variance() {
    ///
    /// (Ex2 - (Ex*Ex)/n) / (n-1)
    ///
    double variance = (this->sum_of_squares - double(this->sum*this->sum)/this->n) / (this->n-1);

    return variance;
}

#endif //RUNLENGTH_ANALYSIS_ITERATIVESUMMARYSTATS_HPP
