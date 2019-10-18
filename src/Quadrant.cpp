
#include "Quadrant.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>

using std::ofstream;
using std::make_shared;
using std::runtime_error;
using std::cout;
using std::to_string;


QuadCoordinate::QuadCoordinate()=default;


QuadCoordinate::QuadCoordinate(float x, float y){
    this->x = x;
    this->y = y;
}


BoundingBox::BoundingBox()=default;


BoundingBox::BoundingBox(QuadCoordinate center, float half_size){
    this->center = center;
    this->half_size = half_size;
}
