//
// Created by jkr on 12/26/24.
//
#include "IndexForce.h"

dvec3 IndexForce::applyForce(Particle &p, const double sim_time) const {
  if (sim_time < time_) {
    return force_values_;
  }
  return {0.0, 0.0, 0.0};
}