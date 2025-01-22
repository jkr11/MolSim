//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

// Disgusting O(n * k * l) search
void IndexForce::applyForce(Particle &p, double sim_time) const {
  if (time < sim_time) {
    p.addF(force_values);
  }
}