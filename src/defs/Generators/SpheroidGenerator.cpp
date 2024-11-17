//
// Created by jkr on 11/17/24.
//

#include "SpheroidGenerator.h"

#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/SpdWrapper.h"
SpheroidGenerator::SpheroidGenerator(const dvec3& origin, const int radius,
                                     const double h, const double m,
                                     const dvec3& initialVelocity,
                                     const double epsilon, const double sigma,
                                     const int type, const bool threeD)
    : origin(origin),
      radius(radius),
      h(h),
      m(m),
      initialVelocity(initialVelocity),
      epsilon(epsilon),
      sigma(sigma),
      type(type),
      threeD(threeD) {
  SpdWrapper::get()->info("CuboidGenerator created with parameters:");
  SpdWrapper::get()->info("origin: ({}, {}, {})", origin[0], origin[1],
                          origin[2]);
  SpdWrapper::get()->info("radius: {}", radius);
  SpdWrapper::get()->info("h: {}", h);
  SpdWrapper::get()->info("m: {}", m);
  SpdWrapper::get()->info("initialVelocity: ({}, {}, {})", initialVelocity[0],
                          initialVelocity[1], initialVelocity[2]);
  SpdWrapper::get()->info("mv: {}", mv);
  SpdWrapper::get()->info("epsilon: {}", epsilon);
  SpdWrapper::get()->info("sigma: {}", sigma);
  SpdWrapper::get()->info("type: {}", type);
}

void SpheroidGenerator::generate(std::vector<Particle>& particles) {
  const int size = static_cast<int>(
      0.75 * M_PI *
      std::pow(radius, 3));  // 1 sphere is always contained in 2 sphere so this
                             // should be fine actually
  particles.reserve(size);
  for (int i = -radius; i <= radius; i++) {
    for (int j = -radius; j <= radius; j++) {
      for (int k = -radius; k <= radius; k++) {
        if (threeD && k == 0) {
          const double spaceRadius = radius * h;
          dvec3 point = {i * h, j * h, k * h};
          if (const double dist = ArrayUtils::L2Norm((1 / spaceRadius) * point);
              dist <= 1.0) {
            dvec3 position = origin + point;
            dvec3 V =
                initialVelocity + maxwellBoltzmannDistributedVelocity(mv, 3);
            particles.emplace_back(position, V, m, epsilon, sigma, type);
          }
        }
      }
    }
  }
}