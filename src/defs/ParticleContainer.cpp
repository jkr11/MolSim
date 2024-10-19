//
// Created by jkr on 10/18/24.
//
#include "ParticleContainer.h"

#include <iostream>
#include <ostream>

ParticleContainer::ParticleContainer() = default;

ParticleContainer::~ParticleContainer() = default;


ParticleContainer::ParticleContainer(const std::list<Particle>& particles) {
  for (const auto & particle : particles) {
    this->particles.push_back(particle);
  }
}

ParticleContainer::ParticleContainer(const std::vector<Particle>& particles) {
  for (const auto & particle : particles) {
    this->particles.push_back(particle);
  }
}



