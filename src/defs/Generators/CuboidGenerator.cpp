//
// Created by jkr on 10/31/24.
//
#include "CuboidGenerator.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/SpdWrapper.h"

CuboidGenerator::CuboidGenerator(const dvec3 &corner,
                                 const std::array<int, 3> &dimensions,
                                 const double h, const double m,
                                 const std::array<double, 3> &initial_velocity,
                                 const double mv, const double epsilon,
                                 const double sigma, const int type,
                                 const bool two_d)
    : corner_(corner),
      dimensions_(dimensions),
      h_(h),
      m_(m),
      initial_velocity_(initial_velocity),
      mv_(mv),
      epsilon_(epsilon),
      sigma_(sigma),
      type_(type),
      two_d_(two_d) {
  DEBUG_PRINT_FMT("CuboidGenerator of dim {} created with parameters:",
                  two_d ? 2 : 3);
  DEBUG_PRINT_FMT("corner: ({}, {}, {})", corner[0], corner[1], corner[2]);
  DEBUG_PRINT_FMT("dimensions: ({}, {}, {})", dimensions[0], dimensions[1],
                  dimensions[2]);
  DEBUG_PRINT_FMT("h: {}", h);
  DEBUG_PRINT_FMT("m: {}", m);
  DEBUG_PRINT_FMT("initialVelocity: ({}, {}, {})", initial_velocity_[0],
                  initial_velocity_[1], initial_velocity_[2]);
  DEBUG_PRINT_FMT("mv: {}", mv);
  DEBUG_PRINT_FMT("epsilon: {}", epsilon);
  DEBUG_PRINT_FMT("sigma: {}", sigma);
  DEBUG_PRINT_FMT("type: {}", type);
}

void CuboidGenerator::generate(std::vector<Particle> &particles) {
  const int size = dimensions_[0] * dimensions_[1] * dimensions_[2];
  const std::size_t offset = particles.size();
  particles.reserve(offset + size);
  DEBUG_PRINT("reserved: " + std::to_string(size) + "particles");

  for (int i = 0; i < dimensions_[0]; i++) {
    for (int j = 0; j < dimensions_[1]; j++) {
      for (int k = 0; k < dimensions_[2]; k++) {
        dvec3 position = {corner_[0] + i * h_, corner_[1] + j * h_,
                          corner_[2] + k * h_};

        std::array<double, 3> v =
            initial_velocity_ +
            maxwellBoltzmannDistributedVelocity(mv_, two_d_ ? 2 : 3);
        particles.emplace_back(position, v, m_, epsilon_, sigma_, type_);
      }
    }
  }
  DEBUG_PRINT("particles: " + std::to_string(particles.size()));
}
