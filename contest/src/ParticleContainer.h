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
  int particleCount = 0;
  ParticleTypeInfo particleTypeInfo;
  LinkedCellsConfig linkedCellsConfig;

  ivec3 cell_count;
  dvec3 cell_dim;
  std::vector<int> cellDivisions;

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

  int getNewId();


};