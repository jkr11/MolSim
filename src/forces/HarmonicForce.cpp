//
// Created by jkr on 12/23/24.
//

#include "HarmonicForce.h"

#include "utils/ArrayUtils.h"

dvec3 HarmonicForce::applyForce(const Particle& p) const {
  dvec3 force_acc = {0.0, 0.0, 0.0};
  for (const auto& [fst, snd] : p.getNeighbours()) {
    dvec3 rv = snd->getX() - p.getX();
    const double r = ArrayUtils::L2Norm(rv);
    force_acc = force_acc + k * (r - (fst ? sr_0 : r_0)) * ((1.0 / r) * rv);
  }
  return force_acc;
}
