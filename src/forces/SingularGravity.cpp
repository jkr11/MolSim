//
// Created by jkr on 11/24/24.
//
#include "SingularGravity.h"

SingularGravity::SingularGravity(const double g) : g_(g){};

dvec3 SingularGravity::applyForce(const Particle& p) const {
  return {0.0, g_ * p.getM(), 0.0};
}