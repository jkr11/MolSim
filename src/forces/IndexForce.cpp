//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

#include "utils/SpdWrapper.h"

// TODO: time
// Disgusting O(n * k * l) search
void IndexForce::applyForce(Particle &p, double sim_time) const {
  p.addF(force_values);
}