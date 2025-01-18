//
// Created by jkr on 12/26/24.
//

#include "IndexForce.h"

#include "utils/SpdWrapper.h"

void IndexForce::applyForce(std::vector<Particle>& particles) {
  SpdWrapper::get()->info("Applying force...");
  for (const auto idx : indeces) {
    particles[idx].addF(force_values);
  }
}
