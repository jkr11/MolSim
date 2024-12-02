//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "Force.h"
/**
 * @brief class exposing calculation of the lennard jones force
 */
class LennardJones final : public Force {
 public:
  LennardJones() = default;
  /**
   * @brief Calculates the lennard-jones force between two particles.
   * @param p1 particle 1
   * @param p2 particle 2
   * @return the singed vector force between particles p1 and p2
   */
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;

  /**
   * @brief calculates the force of the ghost particle
   * @param p Particle to calculate Force for
   * @param distance the distance to the boundary
   * @return the force in just one dimension
   */
  static double simpleForce(const Particle& p, double distance);
};
#endif  // LENNARDJONES_H
