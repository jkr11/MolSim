//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef PARTICLEGENERATOR_H
#define PARTICLEGENERATOR_H
#include "ParticleContainer.h"

// TODO
class ParticleGenerator {
public:
  virtual ~ParticleGenerator() = default;

  /**
   * @brief spawns particles in a given shape
   * @param particles vector of particles
   */
  virtual void generate(std::vector<Particle> &particles) = 0;
};

#endif  // PARTICLEGENERATOR_H
