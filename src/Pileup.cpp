#include "Pileup.hpp"

Pileup::Pileup(size_t n_channels) {
    this->n_channels = n_channels;
    this->pileup.resize(n_channels);
}
