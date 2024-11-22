//
// Created by jkr on 10/18/24.
//
#pragma once
#include "Force.h"
#include "defs/Particle.h"
#include "utils/ArrayUtils.h"

/**
 * @brief Gravitational force F_{1,2}
 */
class Gravity final : public Force {
 public:
  /**
   *
   */
  Gravity() = default;

  /**
   * Compute the gravitational force between two particles p1 and p2
   * @param p1 Particle1
   * @param p2 Particle2
   * @return Force-vector
   */
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};
