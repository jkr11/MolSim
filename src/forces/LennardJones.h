//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "BidirectionalForce.h"
/**
 * @brief class exposing calculation of the lennard jones force
 */
class LennardJones final : public BidirectionalForce {
 public:
  LennardJones() = default;
  /**
   * @brief Calculates the lennard-jones force between two particles.
   * @param p1 particle 1
   * @param p2 particle 2
   * @return the singed vector force between particles p1 and p2
   */
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};
#endif  // LENNARDJONES_H
