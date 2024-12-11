//
// Created by jkr on 10/18/24.
//
#include "Gravity.h"

dvec3 Gravity::directionalForce(Particle& p1, Particle& p2, const double& r_sq,
                       const dvec3& rv) const {
  // const dvec3 rv = p2.getX() - p1.getX();
  //  const double dist = ArrayUtils::L2Norm(r);
  const double F = (p1.getM() * p2.getM()) / r_sq;
  return F / r_sq * (rv);
}

double Gravity::simpleForce(Particle& p, double distance) { return 0; }
