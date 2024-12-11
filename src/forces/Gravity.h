//
// Created by jkr on 10/18/24.
//
#pragma once
#include "InteractiveForce.h"
#include "defs/Particle.h"
#include "utils/ArrayUtils.h"

/**
 * @brief Gravitational force F_{1,2}
 */
class Gravity final : public InteractiveForce {
 public:
  Gravity() = default;

  /**
   * Compute the gravitational force between two particles p1 and p2
   * @param p1 Particle1
   * @param p2 Particle2
   * @return Force-vector
   */
  dvec3 directionalForce(Particle& p1, Particle& p2,
                         const double r) const override;

  /**
   * @brief calculates the force of the ghost particle, not implemented for this
   * Force
   * @param p Particle to calculate Force for
   * @param distance the distance to the boundary
   * @return the force in just one dimension
   */
  static double simpleForce(Particle& p, double distance);
};
