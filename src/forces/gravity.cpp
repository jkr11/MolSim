#include "gravity.h"

// This is the force in direction 1 -> 2, so F_{12}
dvec3 Gravity::directionalForce(Particle& p1, Particle& p2) const {
  const dvec3 r = p2.getX() - p1.getX();
  const double dist = ArrayUtils::L2Norm(r);
  if (dist == 0) {
    return {0, 0, 0};
  }
  const double F = (p1.getM() * p2.getM()) / std::pow(dist, 2);
  return F / dist * (p2.getX() - p1.getX());
}
