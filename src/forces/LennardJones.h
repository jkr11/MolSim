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
   * @return the singed vector force between particles p1 and p2
   */
  dvec3 inline directionalForce(Particle& p1, Particle& p2,
                                const double r_sq) const override {
    const dvec3 rv = p2.getX() - p1.getX();
    // const double r = ArrayUtils::L2Norm(rv);
    const double sigma_sq = std::pow((p1.getSigma() + p2.getSigma()) / 2, 2);
    const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    const double sr_sq = sigma_sq / r_sq;
    const double sr6 = std::pow(sr_sq, 3);
    const double sr12 = std::pow(sr6, 2);
    const double force_magnitude = 24 * epsilon * (sr6 - 2 * sr12) / r_sq;
    return force_magnitude * rv;
  }

  /**
   * @brief calculates the force of the ghost particle
   * @param p Particle to calculate Force for
   * @param r_sq the distance to the boundary
   * @return the force in just one dimension
   */
  static double simpleForce(const Particle& p, double r_sq);
};
#endif  // LENNARDJONES_H
