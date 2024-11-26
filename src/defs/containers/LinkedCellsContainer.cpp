//
// created by mcarn on 11/15/24
//
#include "LinkedCellsContainer.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "debug/debug_print.h"
#include "defs/Particle.h"
#include "utils/SpdWrapper.h"

LinkedCellsContainer::LinkedCellsContainer(const ivec3 &domain,
                                           const double cutoff) {
  DEBUG_PRINT("LinkedCellsContainer instantiated");
  SpdWrapper::get()->info("domain size: ({}, {}, {})", domain[0], domain[1],
                          domain[2]);

  cells = {};
  this->cutoff = cutoff;

  cellCount = {std::max(static_cast<int>(std::floor(domain[0] / cutoff)), 1),
               std::max(static_cast<int>(std::floor(domain[1] / cutoff)), 1),
               std::max(static_cast<int>(std::floor(domain[2] / cutoff)), 1)};

  cellDim = {domain[0] / cellCount[0], domain[1] / cellCount[1],
             domain[2] / cellCount[2]};
  // add 2 for halo

  cellCount = {cellCount[0] + 2, cellCount[1] + 2, cellCount[2] + 2};

  cells.resize(cellCount[0] * cellCount[1] * cellCount[2]);

  // TODO: pretty
  SpdWrapper::get()->info("cell dim: {}, {}, {}; cell count: {}, {}, {}",
                          cellDim[0], cellDim[1], cellDim[2], cellCount[0],
                          cellCount[1], cellCount[2]);
}

void LinkedCellsContainer::addParticle(const Particle &p) {
  // SpdWrapper::get()->info("addParticle");
  const std::size_t index = dvec3ToCellIndex(p.getX());
  cells[index].emplace_back(p);

  DEBUG_PRINT_FMT("Added particle with coords ({}, {}, {}) into cell index: {}",
                  p.getX()[0], p.getX()[1], p.getX()[2], index)
}

void LinkedCellsContainer::addParticles(
    const std::vector<Particle> &particles) {
  for (const Particle &p : particles) {
    addParticle(p);
  }
}

void LinkedCellsContainer::removeParticle(const Particle &p) {
  const std::size_t index = dvec3ToCellIndex(p.getX());
  std::vector<Particle> &particles = cells[index];

  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [&p](const Particle &q) { return p == q; }),
                  particles.end());

  DEBUG_PRINT_FMT(
      "Removed particle with coords ({}, {}, {}) from cell index: {}",
      p.getX()[0], p.getX()[1], p.getX()[2], index)
}

std::vector<Particle *> LinkedCellsContainer::getParticles() {
  std::vector<Particle *> res;
  singleIterator([&res](Particle &p) { res.push_back(&p); });

  return res;
}

[[nodiscard]] std::size_t LinkedCellsContainer::size() const {
  std::size_t count = 0;
  for (auto &c : cells) {
    count += c.size();
  }
  return count;
}

void LinkedCellsContainer::imposeInvariant() {
  for (std::size_t index = 0; index < cells.size(); index++) {
    for (auto it = cells[index].begin(); it < cells[index].end();) {
      const std::size_t shouldBeIndex = dvec3ToCellIndex(it->getX());
      if (shouldBeIndex == index) {
        ++it;
        continue;
      }

      cells[shouldBeIndex].push_back(*it);
      it = cells[index].erase(it);
    }
  }
}

void LinkedCellsContainer::singleIterator(
    const std::function<void(Particle &)> &f) {
  for (auto &c : cells) {
    for (auto &p : c) {
      f(p);
    }
  }
}

void LinkedCellsContainer::pairIterator(
    const std::function<void(Particle &, Particle &)> &f) {
  // - as x, y, z are all increasing offsets point to all neighbours in
  // positive directions
  // - for better cache usage (in flattened version) minimize max distance
  //   between values
  // - direction of traversal: z, y, x
  const std::array<ivec3, 13> offsets = {{
      // 9 x facing
      {{1, -1, -1}},
      {{1, -1, 0}},
      {{1, -1, 1}},
      {{1, 0, -1}},
      {{1, 0, 0}},
      {{1, 0, 1}},
      {{1, 1, -1}},
      {{1, 1, 0}},
      {{1, 1, 1}},
      // 3 y
      {{0, 1, -1}},
      {{0, 1, 0}},
      {{0, 1, 1}},
      // last z
      {{0, 0, 1}},
  }};

  // go over all cell indices
  for (std::size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
    std::vector<Particle> &cellParticles = cells[cellIndex];

    if (cellParticles.empty()) continue;

    ivec3 cellCoordinate = cellIndexToCoord(cellIndex);
    DEBUG_PRINT_FMT("cell index: {}; coord = ({}, {}, {}); halo? = {}",
                    cellIndex, cellCoordinate[0], cellCoordinate[1],
                    cellCoordinate[2], isHalo(cellIndex));

    // iterate over particles inside cell

    for (std::size_t i = 0; i < cellParticles.size(); ++i) {
      for (std::size_t j = i + 1; j < cellParticles.size(); ++j) {
        const dvec3 p = cellParticles[i].getX();
        const dvec3 q = cellParticles[j].getX();
        if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
          continue;
        f(cellParticles[i], cellParticles[j]);
        DEBUG_PRINT_FMT("Intra cell pair: ({}, {})", cellParticles[i].getType(),
                        cellParticles[j].getType());
      }
    }

    // iterate over neighbouring particles
    for (auto &offset : offsets) {
      // compute neighbourIndex and check if it is valid
      const ivec3 neighbourCoord = {cellCoordinate[0] + offset[0],
                                    cellCoordinate[1] + offset[1],
                                    cellCoordinate[2] + offset[2]};

      if (!isValidCellCoordinate(neighbourCoord)) {
        DEBUG_PRINT_FMT("Invalid coord: ({}, {}, {})", neighbourCoord[0],
                        neighbourCoord[1], neighbourCoord[2])
        continue;
      }

      const size_t neighbourIndex = cellCoordToIndex(neighbourCoord);
      DEBUG_PRINT_FMT(
          "Checking cell i={}; c=({}, {}, {}) for pairs (offset = ({}, {}, "
          "{}))",
          neighbourIndex, neighbourCoord[0], neighbourCoord[1],
          neighbourCoord[2], offset[0], offset[1], offset[2]);

      // go over all pairs with neighbour particles
      std::vector<Particle> &neighbourParticles = cells[neighbourIndex];
      if (neighbourParticles.empty()) continue;

      for (auto &cellParticle : cellParticles) {
        for (auto &neighbourParticle : neighbourParticles) {
          auto p = cellParticle.getX();
          auto q = neighbourParticle.getX();

          if (dvec3 d = {p[0] - q[0], p[1] - q[1], p[2] - q[2]};
              d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
            continue;

          f(cellParticle, neighbourParticle);
          DEBUG_PRINT_FMT("Cross cell pair: ({}, {})", cellParticle.getType(),
                          neighbourParticle.getType())
        }
      }
    }
  }
}

void LinkedCellsContainer::boundaryIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells.size(); index++) {
    if (!isBoundary(index)) continue;

    for (auto &p : cells[index]) {
      f(p);
    }
  }
}

void LinkedCellsContainer::haloIterator(
    const std::function<void(Particle &)> &f) {
  for (std::size_t index = 0; index < cells.size(); index++) {
    if (!isHalo(index)) continue;

    for (auto &p : cells[index]) {
      f(p);
    }
  }
}

inline std::size_t LinkedCellsContainer::dvec3ToCellIndex(
    const dvec3 &position) const {
  const std::array<int, 3> cellCoords = {
      static_cast<int>(std::floor(position[0] / cellDim[0])),
      static_cast<int>(std::floor(position[1] / cellDim[1])),
      static_cast<int>(std::floor(position[2] / cellDim[2]))};

  return cellCoordToIndex(cellCoords);
}

inline std::size_t LinkedCellsContainer::cellCoordToIndex(
    const ivec3 position) const {
  return (position[0] + 1) * (cellCount[1] * cellCount[2]) +
         (position[1] + 1) * cellCount[2] + (position[2] + 1);
}

inline ivec3 LinkedCellsContainer::cellIndexToCoord(
    std::size_t cellIndex) const {
  const int x = static_cast<int>(cellIndex / (cellCount[1] * cellCount[2]));
  cellIndex = cellIndex - (x * cellCount[1] * cellCount[2]);

  const int y = static_cast<int>(cellIndex / cellCount[2]);
  const int z = static_cast<int>(cellIndex - (y * cellCount[2]));

  return {x - 1, y - 1, z - 1};
}

inline bool LinkedCellsContainer::isValidCellCoordinate(
    const ivec3 coordinate) const {
  return (-1 <= coordinate[0] && coordinate[0] <= (cellCount[0] - 2)) &&
         (-1 <= coordinate[1] && coordinate[1] <= (cellCount[1] - 2)) &&
         (-1 <= coordinate[2] && coordinate[2] <= (cellCount[2] - 2));
}

inline bool LinkedCellsContainer::isHalo(const ivec3 cellCoord) const {
  return cellCoord[0] == -1 || cellCoord[1] == -1 || cellCoord[2] == -1 ||
         cellCoord[0] == (cellCount[0] - 2) ||
         cellCoord[1] == (cellCount[1] - 2) ||
         cellCoord[2] == (cellCount[2] - 2);
}

inline bool LinkedCellsContainer::isHalo(const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isHalo(cellCoord);
}

inline bool LinkedCellsContainer::isBoundary(const ivec3 cellCoord) const {
  return (cellCoord[0] == 0 || cellCoord[1] == 0 || cellCoord[2] == 0 ||
          cellCoord[0] == (cellCount[0] - 3) ||
          cellCoord[1] == (cellCount[1] - 3) ||
          cellCoord[2] == (cellCount[2] - 3)) &&
         !isHalo(cellCoord);
}

inline bool LinkedCellsContainer::isBoundary(
    const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isBoundary(cellCoord);
}