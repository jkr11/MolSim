//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

#include "debug/debug_print.h"

// TODO: time
// Disgusting O(n * k * l) search
dvec3 IndexForce::applyForce(Particle &p, double sim_time) const {
  InfoVec("Force vals", force_values);
  if (sim_time < time) {
    return force_values;
  }

  return {0.0, 0.0, 0.0};
}