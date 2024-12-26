//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

void IndexForce::applyForce(std::vector<Particle>& particles) {
  for (const auto idx : indeces) {
    particles[idx].addF(force_values);
  }
}
