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
                                 const std::array<double, 3> &initialVelocity,
                                 const double mv, const double epsilon,
                                 const double sigma, const int type,
                                 const bool twoD)
    : corner(corner),
      dimensions(dimensions),
      h(h),
      m(m),
      initialVelocity(initialVelocity),
      mv(mv),
      epsilon(epsilon),
      sigma(sigma),
      type(type),
      twoD(twoD) {
  SpdWrapper::get()->info("CuboidGenerator of dim {} created with parameters:",
                          twoD ? 2 : 3);
  SpdWrapper::get()->info("corner: ({}, {}, {})", corner[0], corner[1],
                          corner[2]);
  SpdWrapper::get()->info("dimensions: ({}, {}, {})", dimensions[0],
                          dimensions[1], dimensions[2]);
  SpdWrapper::get()->info("h: {}", h);
  SpdWrapper::get()->info("m: {}", m);
  SpdWrapper::get()->info("initialVelocity: ({}, {}, {})", initialVelocity[0],
                          initialVelocity[1], initialVelocity[2]);
  SpdWrapper::get()->info("mv: {}", mv);
  SpdWrapper::get()->info("epsilon: {}", epsilon);
  SpdWrapper::get()->info("sigma: {}", sigma);
  SpdWrapper::get()->info("type: {}", type);
}

void CuboidGenerator::generate(std::vector<Particle> &particles) {
  const int size = (dimensions[0] * dimensions[1] * dimensions[2]);
  particles.reserve(size);
  DEBUG_PRINT("reserved: " + std::to_string(size) + "particles");

  for (int i = 0; i < dimensions[0]; i++) {
    for (int j = 0; j < dimensions[1]; j++) {
      for (int k = 0; k < dimensions[2]; k++) {
        dvec3 position = {corner[0] + i * h, corner[1] + j * h,
                          corner[2] + k * h};

        std::array<double, 3> V =
            initialVelocity +
            maxwellBoltzmannDistributedVelocity(mv, twoD ? 2 : 3);
        particles.emplace_back(position, V, m, epsilon, sigma, type);
      }
    }
  }
  DEBUG_PRINT("particles: " + particles.size());
}
