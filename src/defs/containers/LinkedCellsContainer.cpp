/*#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "../../utils/SpdWrapper.h"
#include "../Particle.h"
#include "../containers/ParticleContainer.h"

class LinkedCellsContainer : public ParticleContainer {
 private:
  /**
   * x is left - right
   * y is up - down
   * z is back - front
   *
   * position = x * (cellsY * cellsZ) + y * (cellsZ) + z
   *
  std::vector<std::vector<Particle>> particles;
  std::array<int, 3> cellCount;
  dvec3 cellDim;

 public:
  //TODO this is just bad
  LinkedCellsContainer() : ParticleContainer() {

  };
  LinkedCellsContainer(const std::vector<Particle> &particles) : ParticleContainer(particles) {};

  LinkedCellsContainer(dvec3 domain, double cutoff /*boundary types*) {
    cellCount = {std::max(static_cast<int>(std::floor(domain[0] / cutoff)), 1),
                 std::max(static_cast<int>(std::floor(domain[1] / cutoff)), 1),
                 std::max(static_cast<int>(std::floor(domain[2] / cutoff)), 1)};

    cellDim = {domain[0] / cellCount[0], domain[1] / cellCount[1],
               domain[2] / cellCount[2]};

    // add 2 for halo
    cellCount = {cellCount[0] + 2, cellCount[1] + 2, cellCount[2] + 2};

    // SpdWrapper::get()->info()
  }

  ~LinkedCellsContainer() = default;

  void addParticle(const Particle &p) {
    int index = positionToCellIndex(p.getX());
    particles[index].emplace_back(p);
  }

  void removeParticle(const Particle &p) {
    //TODO
  }

  std::vector<Particle> getParticles() {
    //TODO
  }

  [[nodiscard]] std::size_t size() const {
    std::size_t count = 0;
    for (auto &v : particles) {
      count += v.size();
    }
    return count;
  }

  void singleIterator(std::function<void(Particle &)> f) {
    for (auto &v : particles) {
      for (auto &p : v) {
        f(p);
      }
    }
  }

  void pairIterator(std::function<void(Particle &, Particle &)> f) {
    //TODO
  }

  void boundaryIterator(std::function<void(Particle&)> f) {
    //TODO
  }

  void haloIterator(std::function<void(Particle&)> f) {
    //TODO
  }

  int positionToCellIndex(dvec3 position) {
    // add + 1 to each coordinate to skip halo
    int cellX = static_cast<int>(std::floor(position[0] / cellDim[0])) + 1;
    int cellY = static_cast<int>(std::floor(position[1] / cellDim[1])) + 1;
    int cellZ = static_cast<int>(std::floor(position[2] / cellDim[2])) + 1;

    return cellX * (cellCount[1] * cellCount[2]) + cellY * (cellCount[2]) +
           cellZ;
  }
};*/