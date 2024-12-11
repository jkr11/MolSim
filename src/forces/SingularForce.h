//
// Created by jkr on 11/24/24.
//
#ifndef SINGULARFORCE_H
#define SINGULARFORCE_H
#pragma once
#include "defs/Particle.h"
/**
 * @brief instance of forces that are global and not an interaction force
 */
class SingularForce {
 public:
  SingularForce() = default;
  virtual ~SingularForce() = default;

  /**
   * @brief calculates the singular force exerted on the particle
   * @param p the particle to be operated on
   * @return the force vector on the particle p
   */
  [[nodiscard]] inline virtual dvec3 applyForce(const Particle& p) const = 0;
};
#endif  // SINGULARFORCE_H
