//
// Created by jkr on 11/17/24.
//

#include "SpheroidGenerator.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/SpdWrapper.h"

SpheroidGenerator::SpheroidGenerator(const dvec3& origin, const int radius,
                                     const double h, const double m,
                                     const dvec3& initial_velocity,
                                     const double epsilon, const double sigma,
                                     const int type, const double mv,
                                     const bool two_d)
    : origin(origin),
      radius(radius -
             1),  // needs to be minus one because we consider the origin as one
      h(h),
      m(m),
      initial_velocity(initial_velocity),
      epsilon(epsilon),
      sigma(sigma),
      type(type),
      mv(mv),
      two_d(two_d) {
  DEBUG_PRINT_FMT("SpheroidGenerator of dim {} created with parameters:",
                  two_d ? 2 : 3);
  DEBUG_PRINT_FMT("origin: ({}, {}, {})", origin[0], origin[1], origin[2]);
  DEBUG_PRINT_FMT("radius: {}", radius + 1);
  DEBUG_PRINT_FMT("h: {}", h);
  DEBUG_PRINT_FMT("m: {}", m);
  DEBUG_PRINT_FMT("initialVelocity: ({}, {}, {})", initial_velocity[0],
                  initial_velocity[1], initial_velocity[2]);
  DEBUG_PRINT_FMT("mv: {}", mv);
  DEBUG_PRINT_FMT("epsilon: {}", epsilon);
  DEBUG_PRINT_FMT("sigma: {}", sigma);
  DEBUG_PRINT_FMT("type: {}", type);
}

void SpheroidGenerator::generate(std::vector<Particle>& particles) {
  int size{};
  // Approximate the size by area of disk and volume of sphere
  // There are between (r-h)/h, (r+h)/h particles on one bisector, so we
  // overestimate with (r+h)
  if (two_d) {
    size = static_cast<int>(M_PI * std::pow((radius + h) / h, 2));
  } else {
    size = static_cast<int>((4.0 / 3.0) * M_PI * std::pow((radius + h) / h, 3));
  }
  particles.reserve(size);
  DEBUG_PRINT("Reserved " + std::to_string(size));
  for (int i = -radius; i <= radius; i++) {
    for (int j = -radius; j <= radius; j++) {
      for (int k = -radius; k <= radius; k++) {
        if (two_d && k != 0) {
          continue;
        }
        const double space_radius = radius * h;
        dvec3 point = {i * h, j * h, k * h};
        if (const double dist = ArrayUtils::L2Norm(1 / space_radius * point);
            dist <= 1.0) {
          dvec3 position = origin + point;
          dvec3 v = initial_velocity +
                    maxwellBoltzmannDistributedVelocity(mv, two_d ? 2 : 3);
          particles.emplace_back(position, v, m, epsilon, sigma, type);
        }
      }
    }
  }
}
