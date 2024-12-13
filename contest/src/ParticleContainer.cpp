//
// Created by mcarn on 12/11/24.
//

#include <iostream>

#include "ParticleContainer.h"

ParticleContainer::ParticleContainer(const ParticleTypeInfo particleTypeInfo, const LinkedCellsConfig linkedCellsConfig, const std::size_t initialCapacity) {
  this->particleTypeInfo = particleTypeInfo;
  this->linkedCellsConfig = linkedCellsConfig;

  cellCount = {std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[0] / linkedCellsConfig.cutoff_radius)), 1),
                std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[1] / linkedCellsConfig.cutoff_radius)), 1),
                std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[2] / linkedCellsConfig.cutoff_radius)), 1)};

  cellDim = {static_cast<double>(linkedCellsConfig.domain[0]) / cellCount[0],
              static_cast<double>(linkedCellsConfig.domain[1]) / cellCount[1],
              static_cast<double>(linkedCellsConfig.domain[2]) / cellCount[2]};

  cellCount = {cellCount[0] + 2, cellCount[1] + 2, cellCount[2] + 2};

  partitionSizes.resize(cellCount[0] * cellCount[1] * cellCount[2]);
  partitionStart.resize(cellCount[0] * cellCount[1] * cellCount[2]);
  swapPointers.resize(cellCount[0] * cellCount[1] * cellCount[2]);

  partitionStart[0] = 0;
  swapPointers[0] = 0;

  std::cout << "Cell count: " << cellCount[0] << ", " << cellCount[1] << ", " << cellCount[2] << std::endl;
  std::cout << "Cell dim: " << cellDim[0] << ", " << cellDim[1] << ", " << cellDim[2] << std::endl;

  //init aligned vectors
  ids.reserve(initialCapacity);
  types.reserve(initialCapacity);

  px.reserve(initialCapacity);
  py.reserve(initialCapacity);
  pz.reserve(initialCapacity);

  vx.reserve(initialCapacity);
  vy.reserve(initialCapacity);
  vz.reserve(initialCapacity);

  fx.reserve(initialCapacity);
  fy.reserve(initialCapacity);
  fz.reserve(initialCapacity);

  ofx.reserve(initialCapacity);
  ofy.reserve(initialCapacity);
  ofz.reserve(initialCapacity);
}

void ParticleContainer::addParticle(int id, int type, dvec3 position, dvec3 velocity, dvec3 force, dvec3 oldForce) {
  this->ids.emplace_back(id);
  this->types.emplace_back(type);

  px.emplace_back(position[0]);
  py.emplace_back(position[1]);
  pz.emplace_back(position[2]);

  vx.emplace_back(velocity[0]);
  vy.emplace_back(velocity[1]);
  vz.emplace_back(velocity[2]);

  fx.emplace_back(force[0]);
  fy.emplace_back(force[1]);
  fz.emplace_back(force[2]);

  ofx.emplace_back(oldForce[0]);
  ofy.emplace_back(oldForce[1]);
  ofz.emplace_back(oldForce[2]);
}

void ParticleContainer::removeParticle(int id) {
  //TODO:
}

void ParticleContainer::imposeInvariant() {
  //TODO: ifdef DEBUG
  /*int placedWrong = 0;
  for (std::size_t i = 0; i < ids.size(); i++) {
    std::size_t correctPartition = dvec3ToIndex({px[i], py[i], pz[i]});
    std::size_t targetStart = partitionStart[correctPartition];
    std::size_t targetEnd = targetStart + partitionSizes[correctPartition];

    if (!(targetStart <= i && i < targetEnd))
      placedWrong++;
  }

  if (placedWrong > 0)
    std::cout << "placed wrong: " << placedWrong << std::endl;
  */
  //end debug

  // first pass: calc cell divisions and pointers
  // reset counts
  for (std::size_t i = 0; i < cellCount[0] * cellCount[1] * cellCount[2]; i++) {
     partitionSizes[i] = 0;
  }

  // count how many particles are in each partition
  for (std::size_t i = 0; i < ids.size(); i++) {
    std::size_t index = dvec3ToIndex({px[i], py[i], pz[i]});
    partitionSizes[index]++;
  }

  // calculate partitions start indices and copy into swap pointers
  for (std::size_t i = 1; i < cellCount[0] * cellCount[1] * cellCount[2]; i++) {
    partitionStart[i] = partitionStart[i - 1] + partitionSizes[i - 1];
    swapPointers[i] = partitionStart[i];
  }

  // second pass: reshuffle to impose invariant
  // - if in correct position: skip
  // - if not:
  //     - calc target currentPartition (which is not already in correct position
  //     - swap current with misplaced in target partition
  //     -

  //TODO: add ifdef debug
  //int swapsPerformed = 0;

  for (std::size_t i = 0; i < ids.size(); ) {
    std::size_t currentPartition = dvec3ToIndex({px[i], py[i], pz[i]});
    std::size_t targetStart = partitionStart[currentPartition];
    std::size_t targetEnd = targetStart + partitionSizes[currentPartition];

    if (targetStart <= i && i < targetEnd) {
      i++;
      continue;
    }

    // Skip already correctly placed elements in target partition
    while (swapPointers[currentPartition] < targetEnd) {
      std::size_t targetIndex = swapPointers[currentPartition];
      std::size_t targetElementPartition = dvec3ToIndex({px[targetIndex], py[targetIndex], pz[targetIndex]});

      if (targetElementPartition != currentPartition) {
        break;
      }

      swapPointers[currentPartition]++;
    }

    // Now, swapPointers[currentPartition] points to the next valid target for swapping
    if (swapPointers[currentPartition] < targetEnd) {
      std::size_t targetIndex = swapPointers[currentPartition];

      // Perform the swap
      std::swap(ids[i], ids[targetIndex]);
      std::swap(types[i], types[targetIndex]);
      std::swap(px[i], px[targetIndex]);
      std::swap(py[i], py[targetIndex]);
      std::swap(pz[i], pz[targetIndex]);
      std::swap(vx[i], vx[targetIndex]);
      std::swap(vy[i], vy[targetIndex]);
      std::swap(vz[i], vz[targetIndex]);
      std::swap(fx[i], fx[targetIndex]);
      std::swap(fy[i], fy[targetIndex]);
      std::swap(fz[i], fz[targetIndex]);
      std::swap(ofx[i], ofx[targetIndex]);
      std::swap(ofy[i], ofy[targetIndex]);
      std::swap(ofz[i], ofz[targetIndex]);

      //TODO: add ifdef debug?
      //swapsPerformed++;
      //std::cout << "Swap performed" << std::endl;

      // Increment swapPointers for the target partition
      swapPointers[currentPartition]++;
      continue;
    }

    // If no valid target was found, just move to the next element
    i++;
  }

  //TODO: add ifdef debug
  /*if (swapsPerformed > 0)
    std::cout << "Swaps performed: " << swapsPerformed << std::endl;

  //TODO: add ifdef debug
  std::size_t last = 0;
  for (std::size_t i = 0; i < ids.size(); i++) {
    std::size_t partition = dvec3ToIndex({px[i], py[i], pz[i]});
    if (last > partition)
      std::cout << "Impose invariant: dvec3ToIndex's of particles not increasing." << std::endl;

    last = partition;
  }
   */
  //end debug

}