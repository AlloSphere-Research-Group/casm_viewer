#ifndef INCLUDE_SLICE_HPP
#define INCLUDE_SLICE_HPP

#include "al/core.hpp"
#include "al/core/math/al_Vec.hpp"
#include "al/core/math/al_Mat.hpp"

#include <vector>
#include <algorithm>

using namespace  al;

struct NeighborInfo
{
  size_t index;
  double distance;
};

al::Vec3f find_slicing_dir(const std::vector<al::Vec3f>& v,
                           size_t num_samples = 500,
                           float sample_radius = 0.2,
                           size_t MAX_NEIGHBORS = 8);

std::vector<al::Vec3f> rotate_to_align_to_z(const std::vector<al::Vec3f>& data,
                                            al::Vec3d layer_dir,
                                            al::Vec3f center_of_rotation);

// p: z-axis aligned data
float find_layer_distance(const std::vector<Vec3f>& p);

float findDistanceNormal(const std::vector<float> &data, al::Vec3f normal);

#endif
