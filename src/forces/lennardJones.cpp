//
// Created by jkr on 10/25/24.
//
#include "lennardJones.h"
#include "../defs/Particle.h"
#include "../utils/ArrayUtils.h"

// TODO: implement this, is this much abstraciton necesssary?
dvec3 LennardJonesForce::directionalForce(Particle& p1, Particle& p2) const {
  const dvec3 rv = p2.getX() - p1.getX();
  double r = ArrayUtils::L2Norm(r);
  if (r == 0) return {0, 0, 0};
  constexpr double sr = p1.getSigma()/r;
  return {0, 0, 0};
}
