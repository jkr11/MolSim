//
// Created by jkr on 10/31/24.
//
#include "CuboidGenerator.h"

#include <iostream>
#include <ostream>

#include "../utils/ArrayUtils.h"
#include "../utils/MaxwellBoltzmannDistribution.h"
#include "../utils/SpdWrapper.h"

CuboidGenerator::CuboidGenerator(const dvec3 &corner,
                                 const std::array<int, 3> &dimensions,
                                 const double h, const double m,
                                 const std::array<double, 3> &initialVelocity,
                                 const double mv, const double epsilon, const double sigma, const int type)
    : corner(corner),
      dimensions(dimensions),
      h(h),
      m(m),
      initialVelocity(initialVelocity),
      mv(mv),
      epsilon(epsilon),
      sigma(sigma),
      type(type) {}

void CuboidGenerator::generate(ParticleContainer &container) const {
  std::vector<Particle> particles = container.getParticlesReference();
  particles.reserve(
      static_cast<int>(dimensions[0] * dimensions[1] * dimensions[2]));
  int res = dimensions[0] * dimensions[1] * dimensions[2];
  SpdWrapper::get()->debug("reserved: {}", res);
  for (int i = 0; i < dimensions[0]; i++) {
    for (int j = 0; j < dimensions[1]; j++) {
      for (int k = 0; k < dimensions[2]; k++) {
        dvec3 position = {corner[0] + i * h, corner[1] + j * h,
                          corner[2] + k * h};

        std::array<double, 3> V =
            initialVelocity + maxwellBoltzmannDistributedVelocity(mv, 3);
        Particle p(position, V, m, epsilon, sigma, type);
        container.addParticle(p);
        // particles.emplace_back(position, initialVelocity, m, 1.0, 1.0, type);
      }
    }
  }
  // xcontainer.setParticles(particles);
  SpdWrapper::get()->debug("particles: {}", particles.size());
}