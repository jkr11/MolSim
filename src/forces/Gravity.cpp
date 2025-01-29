//
// Created by jkr on 10/18/24.
//
#include "Gravity.h"

dvec3 Gravity::directionalForce(Particle& p1, Particle& p2) const {
  const dvec3 r = p2.getX() - p1.getX();
  const double dist = ArrayUtils::L2Norm(r);
  const double f = (p1.getM() * p2.getM()) / std::pow(dist, 2);
  return f / dist * (p2.getX() - p1.getX());
}

double Gravity::simpleForce(Particle& p, double distance) { return 0; }
