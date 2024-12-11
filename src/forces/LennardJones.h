//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "InteractiveForce.h"
#include "utils/ArrayUtils.h"  // do not remove this even if clion copes about it
/**
 * @brief class exposing calculation of the lennard jones force
 */
class LennardJones final : public InteractiveForce {
 public:
  LennardJones() = default;
  /**
   * @brief Calculates the lennard-jones force between two particles.
   * @param p1 particle 1
   * @param p2 particle 2
   * @param rv
   * @return the singed vector force between particles p1 and p2
   */
  dvec3 directionalForce(Particle& p1, Particle& p2, const double& r_sq,
                         const dvec3& rv) const override;

  /**
   * @brief calculates the force of the ghost particle
   * @param p Particle to calculate Force for
   * @param r_sq the distance to the boundary
   * @return the force in just one dimension
   */
  static double simpleForce(const Particle& p, double r_sq);
};
#endif  // LENNARDJONES_H
