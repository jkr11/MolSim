//
// Created by mcarn on 11/19/24.
//
#include "FlattenedLinkedCellsContainer.h"

#include "defs/Particle.h"
#include "utils/SpdWrapper.h"

FlattenedLinkedCellsContainer::FlattenedLinkedCellsContainer(
    const ivec3 &domain, const double cutoff) {
  SpdWrapper::get()->info("domain size: ({}, {}, {})", domain[0], domain[1],
                          domain[2]);

  particles = {};
  this->cutoff = cutoff;

  cellCount = {std::max(static_cast<int>(std::floor(domain[0] / cutoff)), 1),
               std::max(static_cast<int>(std::floor(domain[1] / cutoff)), 1),
               std::max(static_cast<int>(std::floor(domain[2] / cutoff)), 1)};

  cellDim = {domain[0] / cellCount[0], domain[1] / cellCount[1],
             domain[2] / cellCount[2]};

  // add 2 for halo
  cellCount = {cellCount[0] + 2, cellCount[1] + 2, cellCount[2] + 2};

  // TODO: pretty
  SpdWrapper::get()->info("cell dim: {}, {}, {}; cell count: {}, {}, {}",
                          cellDim[0], cellDim[1], cellDim[2], cellCount[0],
                          cellCount[1], cellCount[2]);
}

void FlattenedLinkedCellsContainer::addParticle(const Particle &p) {
  const std::size_t index = dvec3ToCellIndex(p.getX());

  // TODO: needs to resize etc, maybe make a addlist that after imposeInvariant
  // adds all automatically
}

void FlattenedLinkedCellsContainer::removeParticle(const Particle &p) {
  // TODO:
}

std::vector<Particle *> FlattenedLinkedCellsContainer::getParticles() {
  // TODO
}

[[nodiscard]] std::size_t FlattenedLinkedCellsContainer::size() const {
  return particles.size();  // TODO: will change if 2 vectors are used to swap
}

void FlattenedLinkedCellsContainer::imposeInvariant() {
  // TODO
}

void FlattenedLinkedCellsContainer::singleIterator(
    const std::function<void(Particle &)> &f) {
  for (auto &p : particles) {
    f(p);
  }
}

void FlattenedLinkedCellsContainer::pairIterator(
    const std::function<void(Particle &, Particle &)> &f) {
  // TODO
}

void FlattenedLinkedCellsContainer::boundaryIterator(
    const std::function<void(Particle &)> &f) {
  // TODO
}

void FlattenedLinkedCellsContainer::haloIterator(
    const std::function<void(Particle &)> &f) {
  // TODO
}

inline std::size_t FlattenedLinkedCellsContainer::dvec3ToCellIndex(
    const dvec3 &position) const {
  const std::array<int, 3> cellCoords = {
      static_cast<int>(std::floor(position[0] / cellDim[0])),
      static_cast<int>(std::floor(position[1] / cellDim[1])),
      static_cast<int>(std::floor(position[2] / cellDim[2]))};

  return cellCoordToIndex(cellCoords);
}

inline std::size_t FlattenedLinkedCellsContainer::cellCoordToIndex(
    const ivec3 position) const {
  // TODO
}

inline ivec3 FlattenedLinkedCellsContainer::cellIndexToCoord(
    std::size_t cellIndex) const {
  // TODO
}

inline bool FlattenedLinkedCellsContainer::isValidCellCoordinate(
    const ivec3 coordinate) const {
  return (-1 <= coordinate[0] && coordinate[0] <= (cellCount[0] - 2)) &&
         (-1 <= coordinate[1] && coordinate[1] <= (cellCount[1] - 2)) &&
         (-1 <= coordinate[2] && coordinate[2] <= (cellCount[2] - 2));
}

inline bool FlattenedLinkedCellsContainer::isHalo(const ivec3 cellCoord) const {
  return cellCoord[0] == -1 || cellCoord[1] == -1 || cellCoord[2] == -1 ||
         cellCoord[0] == (cellCount[0] - 2) ||
         cellCoord[1] == (cellCount[1] - 2) ||
         cellCoord[2] == (cellCount[2] - 2);
}

inline bool FlattenedLinkedCellsContainer::isHalo(
    const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isHalo(cellCoord);
}

inline bool FlattenedLinkedCellsContainer::isBoundary(
    const ivec3 cellCoord) const {
  return (cellCoord[0] == 0 || cellCoord[1] == 0 || cellCoord[2] == 0 ||
          cellCoord[0] == (cellCount[0] - 3) ||
          cellCoord[1] == (cellCount[1] - 3) ||
          cellCoord[2] == (cellCount[2] - 3)) &&
         !isHalo(cellCoord);
}

inline bool FlattenedLinkedCellsContainer::isBoundary(
    const std::size_t cellIndex) const {
  const ivec3 cellCoord = cellIndexToCoord(cellIndex);
  return isBoundary(cellCoord);
}
