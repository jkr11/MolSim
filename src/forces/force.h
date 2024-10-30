//
// Created by jkr on 10/18/24.
//
#pragma once
#include "../defs/Particle.h"

/**
 * @brief Interface for different types of forces
 * @note Force is an interface as we have at least one more "Type" of force
 *       (Lennard-Jones) coming up. Maybe also include self potentials, so split in
 *       (p,q) interactions and self interactions
 */
class Force {
 public:
  /**
   * @brief Create Force object
   * @note Since this is an interface, it's invalid
   */
  Force() = delete;

  /**
   * @brief Virtual destructor for all Force inheritors
   */
  virtual ~Force() = default;

  /**
   * @brief Function to calculate the directional force between two particles
   * @param p1
   * @param p2
   * @return Force-vector
   */
  virtual dvec3 directionalForce(Particle& p1, Particle& p2) const = 0;
};
