//
// Created by mcarn on 12/11/24.
//

#pragma once
#include <vector>
#include <cmath>

#include "../../src/defs/types.h"
#include "Simulation.h"

class ParticleContainer {
public:
  int idCounter = 0;
  ParticleTypeInfo particleTypeInfo;
  LinkedCellsConfig linkedCellsConfig;

  ivec3 cellCount;
  dvec3 cellDim;
  std::vector<std::size_t> partitionSizes;
  std::vector<std::size_t> partitionStart;
  std::vector<std::size_t> swapPointers;

  std::vector<int> ids;
  std::vector<int> types;

  std::vector<double> px, py, pz;
  std::vector<double> vx, vy, vz;
  std::vector<double> fx, fy, fz;
  std::vector<double> ofx, ofy, ofz;

public:
  ParticleContainer(const ParticleTypeInfo particleTypeInfo, const LinkedCellsConfig linkedCellsConfig);
  ~ParticleContainer() = default;

  void addParticle(int id, int type, dvec3 position, dvec3 velocity, dvec3 force, dvec3 oldForce);
  void removeParticle(int id);

  void imposeInvariant();

  //TODO: check if compiler unwraps ivec3 when inlining
  inline std::size_t dvec3ToIndex(const dvec3 &position) const {
    //TODO: add ifdef DEBUG, then check if pos outside of domain

    const std::array<int, 3> cellCoords = {
      static_cast<int>(std::floor(position[0] / cellDim[0])),
      static_cast<int>(std::floor(position[1] / cellDim[1])),
      static_cast<int>(std::floor(position[2] / cellDim[2]))};

    return cellCoordinateToIndex(cellCoords);
  }

  //TODO: check if compiler unwraps ivec3 when inlining
  inline std::size_t cellCoordinateToIndex(const ivec3 &coordinate) const {
    //TODO: add ifdef DEBUG, then check if coordinate outside of domain

    return (coordinate[0] + 1) * (cellCount[1] * cellCount[2]) +
         (coordinate[1] + 1) * (cellCount[2]) + (coordinate[2] + 1);
  };

  //TODO: check if compiler unwraps ivec3 when inlining
  inline ivec3 cellIndexToCoordinate(std::size_t index) const {
    //TODO: add ifdef DEBUG, then check if index outside of cellDomain

    const int x = static_cast<int>(index / (cellCount[1] * cellCount[2]));
    index = index - (x * cellCount[1] * cellCount[2]);

    const int y = static_cast<int>(index / cellCount[2]);
    const int z = static_cast<int>(index - (y * cellCount[2]));

    return {x - 1, y - 1, z - 1};
  };

  inline bool isValidCellCoordinate(
    const ivec3 coordinate) const {
    return (-1 <= coordinate[0] && coordinate[0] <= (cellCount[0] - 2)) &&
           (-1 <= coordinate[1] && coordinate[1] <= (cellCount[1] - 2)) &&
           (-1 <= coordinate[2] && coordinate[2] <= (cellCount[2] - 2));
  }

  inline int getNewId() {
    return ++idCounter;
  }

  //TODO: add abstractions to load a particles position using index, also storing etc. for other attributes ???
  // does this allow us to integrate it more easily into our current main program ???


};