//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef PARTICLEGENERATOR_H
#define PARTICLEGENERATOR_H
#include "ParticleContainer.h"

class ParticleGenerator {
 public:
  virtual ~ParticleGenerator() = default;
  /**
   * @brief spawns particles in a given shape
   * @param container container holding the vector of particles
   */
  virtual void generate(ParticleContainer &container) const = 0;
};

#endif  // PARTICLEGENERATOR_H
