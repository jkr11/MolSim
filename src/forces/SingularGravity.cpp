//
// Created by jkr on 11/24/24.
//
#include "SingularGravity.h"

SingularGravity::SingularGravity(const double g, const int axis)
    : g(g), axis(axis) {};

dvec3 SingularGravity::applyForce(const Particle& p) const {
  dvec3 force{0.0, 0.0, 0.0};
  force[axis] = g * p.getM();
  return force;
}