
#include "Pileup.hpp"

Pileup::Pileup() {
    this->n_channels = 0;
    this->pileup = vector <vector <vector <float> > >();
}


Pileup::Pileup(size_t n_channels, size_t region_size, size_t maximum_depth, vector<float>& default_read_data) {
    this->n_channels = n_channels;
    this->pileup = vector <vector <vector <float> > >(region_size, vector <vector <float> >(maximum_depth, default_read_data));
    this->coverage_per_position.resize(region_size, 0);
}
