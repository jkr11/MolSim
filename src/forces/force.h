//
// Created by jkr on 10/18/24.
//
#pragma once
#include "../defs/Particle.h"

/**
 * @brief Interface for different types of forces
 * @note Force is a class as we have at least one more "Type" of force
 (Lennard-Jones) coming up. Maybe also include self potentials, so split in
 (p,q) interactions and self interactions
 */
class Force {
 public:
  virtual ~Force() = default;

  virtual dvec3 directionalForce(Particle& p1, Particle& p2) const = 0;
};
