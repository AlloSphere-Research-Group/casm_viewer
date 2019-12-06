#include "slice.hpp"
#include "al/math/al_Matrix4.hpp"
#include "al/math/al_Random.hpp"

al::Vec3f find_slicing_dir(const std::vector<al::Vec3f> &v, size_t num_samples,
                           float sample_radius, size_t MAX_NEIGHBORS) {
  al::Vec3f averaged_averaged{0.0, 0.0, 0.0};
  // repeat finding normal direction from random sample points
  for (size_t i = 0; i < num_samples; i += 1) {
    int rand_idx = al::rnd::uniformi(v.size() - 1);
    auto &sample = v[rand_idx];
    // first find all neighbors in range
    std::vector<NeighborInfo> neighbors;
    neighbors.reserve(v.size());
    for (size_t j = 0; j < v.size(); j++) {
      const auto &p = v[j];
      double dist = (p - sample).mag();
      if (0 < dist && dist < sample_radius) {
        neighbors.push_back(NeighborInfo{j, dist});
      }
    }
    if (neighbors.size() == 0) {
      return averaged_averaged;
    }
    // sort the neighbors
    std::sort(neighbors.begin(), neighbors.end(),
              [](const NeighborInfo &a, const NeighborInfo &b) {
                return a.distance < b.distance;
              });
    // cut the number of neighbors
    if (neighbors.size() > MAX_NEIGHBORS) {
      neighbors.resize(MAX_NEIGHBORS);
    }
    // get normals from the sample point and two neighbors
    std::vector<al::Vec3f> normals;
    const auto &c = sample; // the center point
    for (size_t j = 0; j < neighbors.size() - 1; j += 1) {
      const auto &a = v[neighbors[j].index];
      const auto &b = v[neighbors[j + 1].index];
      auto &n = al::cross(a - c, b - c).normalize();
      if (n.z < 0.0f) {
        n *= -1.0f;
      } else if (n.z == 0.0f) {
        if (n.y < 0.0f)
          n *= -1.0f;
      }
      normals.push_back(n);
    }
    // average the normals
    al::Vec3f averaged{0.0, 0.0, 0.0};
    //    std::cout << "Normals ";
    for (const auto &n : normals) {
      //        std::cout << n << " ";
      averaged += n;
    }
    //    std::cout << std::endl;
    averaged /= double(normals.size());
    averaged.normalize();

    //    std::cout << "Averaged: " << averaged << std::endl;
    averaged_averaged += averaged;
  }
  // average the normals from each sample points
  averaged_averaged /= double(num_samples);
  averaged_averaged.normalize();
  return averaged_averaged;
}

std::vector<Vec3f> rotate_to_align_to_z(const std::vector<Vec3f> &data,
                                        Vec3d layer_dir,
                                        Vec3f center_of_rotation) {
  auto v = data; // copy
  for (auto &p : v) {
    al::Mat4f t;
    t.setIdentity();
    // [!] mat mults right to left (M3 * M2 * M1 * v)
    t *= al::Matrix4f::translation(center_of_rotation); // put back
    const auto rot = al::Quatf::getRotationTo(layer_dir, {0.0f, 0.0f, 1.0f});
    al::Mat4f rot_mat;
    rot.toMatrix(rot_mat.elems());
    t *= rot_mat;
    t *=
        al::Matrix4f::translation(-center_of_rotation); // to center of rotation
    p.set(t * al::Vec4f(p, 1)); // need vec4 for multiplying with mat4
  }
  return v;
}

float find_layer_distance(const std::vector<Vec3f> &p) {
  const size_t N = p.size();
  const size_t n = 50; // take top n z-differences
  // if number of layers in data is larger than this
  // result may not be correct

  // get only z values
  std::vector<float> zs;
  zs.resize(N);
  for (size_t i = 0; i < N; i += 1) {
    zs[i] = p[i].z;
  }
  // sort z ascending order
  sort(zs.begin(), zs.end());

  // differences in z: if big, means layer jump
  std::vector<float> dz;
  dz.resize(N - 1);
  for (size_t i = 0; i < N - 1; i += 1) {
    dz[i] = zs[i + 1] - zs[i];
  }
  // sort descending so big jumps come to beginning
  std::sort(dz.begin(), dz.end(), std::greater<float>());

  // in top n z-differences, find maximum decrease and its index.
  // values before this index should be layer distances
  // values after this index should be z-differences between
  // data points in same layer
  size_t max_idx = 0;
  float ddz = 0, ddz_temp = 0;
  for (size_t i = 0; i < n - 1; i += 1) {
    ddz_temp = dz[i] - dz[i + 1];
    if (ddz_temp > ddz) {
      ddz = ddz_temp;
      max_idx = i;
    }
  }

  // average layer distances
  float avg = 0;
  for (size_t i = 0; i < max_idx + 1; i += 1) {
    avg += dz[i];
  }
  avg = avg / (max_idx + 1);

  // above average is minimum distance between layers
  // (dist between max point from layer0 and min point from layer1)
  // needed is average distance between them
  // cluster z value points
  std::vector<float> clustered;
  {
    float z_sum = zs[0];
    int z_cnt = 1;
    for (size_t i = 1; i < N; i += 1) {
      // 0.5 * avg with assumption that spread inside layer is
      // smaller than half layer distance
      if (zs[i] - zs[i - 1] < 0.5 * avg) {
        z_sum += zs[i];
        z_cnt += 1;
      } else {
        z_sum /= z_cnt;
        clustered.push_back(z_sum);
        z_sum = zs[i];
        z_cnt = 1;
      }
    }
    // add last layer
    z_sum /= z_cnt;
    clustered.push_back(z_sum);
  }

  // calc differences between clustered layers and average them
  float result = 0;
  for (size_t i = 0; i < clustered.size() - 1; i += 1) {
    result += clustered[i + 1] - clustered[i];
  }
  result /= (clustered.size() - 1);

  return result;
}

float findDistanceNormal(const std::vector<float> &data, Vec3f normal) {
  const size_t N = data.size() / 4;
  const size_t n = 50; // take top n z-differences
  // if number of layers in data is larger than this
  // result may not be correct

  if (N <= 0) {
    return 1.0;
  }
  normal = normal.normalized();
  // get only z values
  std::vector<float> dist;
  dist.resize(N);
  for (size_t i = 0; i < N; i++) {
    al::Vec3f pos(data[i * 4], data[i * 4 + 1], data[i * 4 + 2]);
    dist[i] = pos.dot(normal);
  }
  // sort z ascending order
  sort(dist.begin(), dist.end());

  // differences in z: if big, means layer jump
  std::vector<float> dz;
  dz.resize(N - 1);
  for (size_t i = 0; i < N - 1; i += 1) {
    dz[i] = dist[i + 1] - dist[i];
  }
  // sort descending so big jumps come to beginning
  std::sort(dz.begin(), dz.end(), std::greater<float>());

  // in top n z-differences, find maximum decrease and its index.
  // values before this index should be layer distances
  // values after this index should be z-differences between
  // data points in same layer
  size_t max_idx = 0;
  float ddz = 0, ddz_temp = 0;
  for (size_t i = 0; i < n - 1; i += 1) {
    ddz_temp = dz[i] - dz[i + 1];
    if (ddz_temp > ddz) {
      ddz = ddz_temp;
      max_idx = i;
    }
  }

  // average layer distances
  float avg = 0;
  for (size_t i = 0; i < max_idx + 1; i += 1) {
    avg += dz[i];
  }
  avg = avg / (max_idx + 1);

  // above average is minimum distance between layers
  // (dist between max point from layer0 and min point from layer1)
  // needed is average distance between them
  // cluster z value points
  std::vector<float> clustered;
  {
    float z_sum = dist[0];
    int z_cnt = 1;
    for (size_t i = 1; i < N; i += 1) {
      // 0.5 * avg with assumption that spread inside layer is
      // smaller than half layer distance
      if (dist[i] - dist[i - 1] < 0.5 * avg) {
        z_sum += dist[i];
        z_cnt += 1;
      } else {
        z_sum /= z_cnt;
        clustered.push_back(z_sum);
        z_sum = dist[i];
        z_cnt = 1;
      }
    }
    // add last layer
    z_sum /= z_cnt;
    clustered.push_back(z_sum);
  }

  // calc differences between clustered layers and average them
  float result = 0;
  for (size_t i = 0; i < clustered.size() - 1; i += 1) {
    result += clustered[i + 1] - clustered[i];
  }
  result /= (clustered.size() - 1);

  return result;
}
