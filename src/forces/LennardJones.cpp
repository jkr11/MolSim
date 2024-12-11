//
// Created by jkr on 10/25/24.
//
#include "LennardJones.h"

#include "utils/ArrayUtils.h"  // do not remove this even if clion copes about it

/*
dvec3 LennardJones::directionalForce(Particle& p1, Particle& p2,
                                     const double r_sq) const {
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
*/

double LennardJones::simpleForce(const Particle& p, const double r_sq) {
  // if (distance <= 0) return 0; shouldnt happen here
  const double sr = std::pow(p.getSigma(),2) / r_sq;
  const double sr6 = std::pow(sr, 3);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * p.getEpsilon() * (sr6 - 2 * sr12) / r_sq;
  return force_magnitude;
}
