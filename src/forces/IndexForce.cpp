//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

#include "debug/debug_print.h"
#include "iostream"

// TODO: time
// Disgusting O(n * k * l) search
dvec3 IndexForce::applyForce(Particle &p, const double sim_time) const {
  // InfoVec("Force vals", force_values_);
  if (sim_time < time_) {
    std::cout << "INdex force at time " << sim_time << "\n";
    return force_values_;
  }

  return {0.0, 0.0, 0.0};
}