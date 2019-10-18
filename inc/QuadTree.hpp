
#ifndef RUNLENGTH_ANALYSIS_QUADTREE_H
#define RUNLENGTH_ANALYSIS_QUADTREE_H

#include <memory>
#include <vector>
#include <array>
#include <string>
#include <experimental/filesystem>
#include "QuadLoss.hpp"
#include "Quadrant.hpp"

using std::vector;
using std::array;
using std::shared_ptr;
using std::string;
using std::experimental::filesystem::path;


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
    vector<QuadCoordinate> points;

    /// Methods ///
    QuadTree();
    QuadTree(BoundingBox boundary);
    bool insert(QuadCoordinate coordinate);
    void subdivide();
    void redistribute_points();
    void write_as_dot(path output_dir="output/", string suffix="0", bool plot=false);
    void append_dot_string(string& edges, string& labels, uint64_t& n);
    void write_bounds(path output_dir="output/", string suffix="0");
    void append_bounds_string(string& bounds);
    void query_range(vector<QuadCoordinate>& results, BoundingBox& bounds);
    uint8_t find_quadrant(QuadCoordinate c);
    virtual shared_ptr<QuadTree> generate_child(BoundingBox bounds);
    virtual void update_loss_from_range(BoundingBox& bounds, QuadLoss& loss_calculator);
};


#endif //RUNLENGTH_ANALYSIS_QUADTREE_H
