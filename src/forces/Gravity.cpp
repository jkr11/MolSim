//
// Created by jkr on 10/18/24.
//
#include "Gravity.h"

dvec3 Gravity::directionalForce(Particle& p1, Particle& p2,
                                const double r) const {
  const dvec3 rv = p2.getX() - p1.getX();
  // const double dist = ArrayUtils::L2Norm(r);
  const double F = (p1.getM() * p2.getM()) / std::pow(r, 2);
  return F / r * (rv);
}

double Gravity::simpleForce(Particle& p, double distance) { return 0; }
