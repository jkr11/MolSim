//
// Created by jkr on 10/18/24.
//
#include "Gravity.h"

dvec3 Gravity::directionalForce(Particle& p1, Particle& p2) const {
  const dvec3 r = p2.getX() - p1.getX();
  const double dist = ArrayUtils::L2Norm(r);

  if (dist == 0) {
    return {0, 0, 0};
  }
  const double F = (p1.getM() * p2.getM()) / std::pow(dist, 2);
  return F / dist * (p2.getX() - p1.getX());
}
