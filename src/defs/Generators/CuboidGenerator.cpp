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
  DEBUG_PRINT_FMT("CuboidGenerator of dim {} created with parameters:",
                  twoD ? 2 : 3);
  DEBUG_PRINT_FMT("corner: ({}, {}, {})", corner[0], corner[1], corner[2]);
  DEBUG_PRINT_FMT("dimensions: ({}, {}, {})", dimensions[0], dimensions[1],
                  dimensions[2]);
  DEBUG_PRINT_FMT("h: {}", h);
  DEBUG_PRINT_FMT("m: {}", m);
  DEBUG_PRINT_FMT("initialVelocity: ({}, {}, {})", initialVelocity[0],
                  initialVelocity[1], initialVelocity[2]);
  DEBUG_PRINT_FMT("mv: {}", mv);
  DEBUG_PRINT_FMT("epsilon: {}", epsilon);
  DEBUG_PRINT_FMT("sigma: {}", sigma);
  DEBUG_PRINT_FMT("type: {}", type);
}

void CuboidGenerator::generate(std::vector<Particle> &particles) {
  const int size = dimensions[0] * dimensions[1] * dimensions[2];
  const std::size_t offset = particles.size();
  particles.reserve(offset + size);
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
#define MEM
#ifdef MEM
  for (int i = 0; i < dimensions[0]; i++) {
    for (int j = 0; j < dimensions[1]; j++) {
      for (int k = 0; k < dimensions[2]; k++) {
        const std::size_t currentIndex =
            i * dimensions[1] * dimensions[2] + j * dimensions[2] + k;

        // Iterate over all neighbors including diagonals
        for (int di = -1; di <= 1; di++) {
          for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
              if (di == 0 && dj == 0 && dk == 0) {
                continue;
              }

              const long ni = i + di;
              const long nj = j + dj;
              const long nk = k + dk;

              if (ni >= 0 && ni < dimensions[0] && nj >= 0 &&
                  nj < dimensions[1] && nk >= 0 && nk < dimensions[2]) {
                const long neighborIndex = ni * dimensions[1] * dimensions[2] +
                                           nj * dimensions[2] + nk;
                const bool isDiagonal = (di != 0) + (dj != 0) + (dk != 0) > 1;
                particles[currentIndex].pushBackNeighbour(
                    isDiagonal, particles[neighborIndex]);
              }
            }
          }
        }
      }
    }
  }
#endif

  DEBUG_PRINT("particles: " + std::to_string(particles.size()));
}
