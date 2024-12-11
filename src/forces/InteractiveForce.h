//
// Created by jkr on 10/18/24.
//
#pragma once
#include "defs/Particle.h"

/**
 * @brief Interface for different types of forces
 * @note Force is an interface as we have at least one more "Type" of force
 *       (Lennard-Jones) coming up. Maybe also include self potentials, so split
 * in (p,q) interactions and self interactions
 */
class InteractiveForce {
 public:
  /**
   * @brief Create Force object
   * @note Since this is an interface, it's invalid
   */
  InteractiveForce() = default;

  /**
   * @brief Virtual destructor for all Force inheritors
   */
  virtual ~InteractiveForce() = default;

  /**
   * @brief Function to calculate the directional force between two particles
   * @param p1
   * @param p2
   * @return Force-vector
   */
  virtual inline dvec3 directionalForce(Particle& p1, Particle& p2,
                                        double r) const = 0;

  /**
   * @brief calculates the force of the ghost particle
   * @param p Particle to calculate Force for
   * @param distance the distance to the boundary
   * @return the force in just one dimension
   */
  static double simpleForce(Particle& p, double distance);
};
