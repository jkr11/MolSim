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
    : origin_(origin),
      radius_(radius - 1),  // needs to be minus one because we consider the
                            // origin as one
      h_(h),
      m_(m),
      initial_velocity_(initial_velocity),
      epsilon_(epsilon),
      sigma_(sigma),
      type_(type),
      mv_(mv),
      two_d_(two_d) {
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
  if (two_d_) {
    size = static_cast<int>(M_PI * std::pow((radius_ + h_) / h_, 2));
  } else {
    size =
        static_cast<int>((4.0 / 3.0) * M_PI * std::pow((radius_ + h_) / h_, 3));
  }
  particles.reserve(size);
  DEBUG_PRINT("Reserved " + std::to_string(size));
  for (int i = -radius_; i <= radius_; i++) {
    for (int j = -radius_; j <= radius_; j++) {
      for (int k = -radius_; k <= radius_; k++) {
        if (two_d_ && k != 0) {
          continue;
        }
        const double space_radius = radius_ * h_;
        dvec3 point = {i * h_, j * h_, k * h_};
        if (const double dist = ArrayUtils::L2Norm(1 / space_radius * point);
            dist <= 1.0) {
          dvec3 position = origin_ + point;
          dvec3 v = initial_velocity_ +
                    maxwellBoltzmannDistributedVelocity(mv_, two_d_ ? 2 : 3);
          particles.emplace_back(position, v, m_, epsilon_, sigma_, type_);
        }
      }
    }
  }
}
