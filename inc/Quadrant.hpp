
#ifndef RUNLENGTH_ANALYSIS_QUADRANT_H
#define RUNLENGTH_ANALYSIS_QUADRANT_H

#include <memory>
#include <vector>
#include <array>
#include <string>
#include <experimental/filesystem>

using std::vector;
using std::array;
using std::shared_ptr;
using std::string;
using std::experimental::filesystem::path;


class QuadCoordinate{
public:
    float x;
    float y;

    QuadCoordinate();
    QuadCoordinate(float x, float y);
};


class BoundingBox{
public:
    QuadCoordinate center;
    float half_size;

    BoundingBox();
    BoundingBox(QuadCoordinate center, float half_size);
    bool contains(QuadCoordinate coordinate);
    bool intersects(BoundingBox boundary);
};


#endif //RUNLENGTH_ANALYSIS_QUADRANT_H
