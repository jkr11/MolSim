//
// Created by mcarn on 12/11/24.
//

#include "ParticleContainer.h"

ParticleContainer::ParticleContainer(const ParticleTypeInfo particleTypeInfo, const LinkedCellsConfig linkedCellsConfig) {
  this->particleTypeInfo = particleTypeInfo;
  this->linkedCellsConfig = linkedCellsConfig;

  cell_count = {std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[0] / linkedCellsConfig.cutoff_radius)), 1),
                std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[1] / linkedCellsConfig.cutoff_radius)), 1),
                std::max(static_cast<int>(std::floor(linkedCellsConfig.domain[2] / linkedCellsConfig.cutoff_radius)), 1)};

  cell_dim = {static_cast<double>(linkedCellsConfig.domain[0]) / cell_count[0],
              static_cast<double>(linkedCellsConfig.domain[1]) / cell_count[1],
              static_cast<double>(linkedCellsConfig.domain[2]) / cell_count[2]};

  cell_count = {cell_count[0] + 2, cell_count[1] + 2, cell_count[2] + 2};


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
  //TODO:
}

int ParticleContainer::getNewId() {
  return ++particleCount;
}