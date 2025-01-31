//
// Created by jkr on 10/31/24.
//
#include "MembraneGenerator.h"

#include <iomanip>
#include <iostream>

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/SpdWrapper.h"

MembraneGenerator::MembraneGenerator(
    const dvec3 &corner, const std::array<int, 3> &dimensions, const double h,
    const double m, const std::array<double, 3> &initial_velocity,
    const double mv, const double epsilon, const double sigma, const int type,
    const bool two_d, const std::vector<ivec3> &indices)
    : corner_(corner),
      dimensions_(dimensions),
      h_(h),
      m_(m),
      initial_velocity_(initial_velocity),
      mv_(mv),
      epsilon_(epsilon),
      sigma_(sigma),
      type_(type),
      two_d_(two_d),
      indices_(indices) {
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

void MembraneGenerator::generate(std::vector<Particle> &particles) {
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

        Particle particle(position, v, m_, epsilon_, sigma_, type_);

        particles.push_back(std::move(particle));
        for (ivec3 vec : indices_) {
          if (ivec3{i, j, k} == vec) {
            SpdWrapper::get()->info("{}, {}, {} matched index {}", i, j, k,
                                    particles.back().getId());
            ids_.push_back(particles.back().getId());
          }
        }
      }
    }
  }

  for (int i = 0; i < dimensions_[0]; i++) {
    for (int j = 0; j < dimensions_[1]; j++) {
      for (int k = 0; k < dimensions_[2]; k++) {
        const std::size_t current_index =
            i * dimensions_[1] * dimensions_[2] + j * dimensions_[2] + k;

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

              if (ni >= 0 && ni < dimensions_[0] && nj >= 0 &&
                  nj < dimensions_[1] && nk >= 0 && nk < dimensions_[2]) {
                const long neighbor_index =
                    ni * dimensions_[1] * dimensions_[2] + nj * dimensions_[2] +
                    nk;
                const bool is_diagonal = (di != 0) + (dj != 0) + (dk != 0) > 1;
                Particle *neighbor_particle = &particles[neighbor_index];
                if (neighbor_particle->getId() == 0) {
                  std::cout
                      << std::hex << std::setw(16) << std::setfill('0')
                      << reinterpret_cast<size_t>(&particles[neighbor_index])
                      << std::dec << std::endl;
                }
                particles[current_index].pushBackNeighbour(
                    is_diagonal, reinterpret_cast<size_t>(neighbor_particle));
              }
            }
          }
        }
      }
    }
  }
  DEBUG_PRINT("particles: " + std::to_string(particles.size()));
}

std::vector<int> MembraneGenerator::getIndices() const { return ids_; }
