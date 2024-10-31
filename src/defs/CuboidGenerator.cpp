//
// Created by jkr on 10/31/24.
//
#include "CuboidGenerator.h"

CuboidGenerator::CuboidGenerator(const dvec3 &corner,
                                 const std::array<int, 3> &dimensions,
                                 const double h, const double m,
                                 const std::array<double, 3> &initialVelocity,
                                 const double temperature, const int type)
    : corner(corner),
      dimensions(dimensions),
      h(h),
      m(m),
      initialVelocity(initialVelocity),
      temp(temperature),
      type(type) {}

void CuboidGenerator::generate(ParticleContainer &container) const {
  std::vector<Particle> particles = container.getParticlesReference();
  particles.reserve(
      static_cast<int>(dimensions[0] * dimensions[1] * dimensions[2]));
  for (int i = 0; i < dimensions[0]; i++) {
    for (int j = 0; j < dimensions[1]; j++) {
      for (int k = 0; k < dimensions[2]; k++) {
        dvec3 position = {corner[0] + i * h, corner[1] + j * h,
                          corner[2] + k * h};
        particles.emplace_back(position, initialVelocity, m, 1.0, 1.0);
      }
    }
  }
}