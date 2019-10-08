
#ifndef RUNLENGTH_ANALYSIS_QUADTREE_H
#define RUNLENGTH_ANALYSIS_QUADTREE_H

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


class QuadTree{
public:
    /// Attributes ///

    static const size_t capacity = 4;

    // Index macros
    static const uint8_t TOP_LEFT = 0;
    static const uint8_t TOP_RIGHT = 1;
    static const uint8_t BOTTOM_LEFT = 2;
    static const uint8_t BOTTOM_RIGHT = 3;
    static const uint8_t NOT_FOUND = -1;

    BoundingBox boundary;

    array<shared_ptr<QuadTree>,4> subtrees;

    /// Methods ///
    QuadTree(BoundingBox boundary);
    bool insert(QuadCoordinate coordinate);
    void subdivide();
    void redistribute_points();
    void find();
    void write_as_dot(path output_dir);
    void append_dot_string(string& edges, string& labels, uint64_t n);

private:
    /// Attributes ///
    vector<QuadCoordinate> points;

    /// Methods ///
    uint8_t find_quadrant(QuadCoordinate c);
};


#endif //RUNLENGTH_ANALYSIS_QUADTREE_H
