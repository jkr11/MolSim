//
// Created by jkr on 12/23/24.
//

#include "HarmonicForce.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

dvec3 HarmonicForce::applyForce(const Particle& p) const {
  dvec3 force_acc = {0.0, 0.0, 0.0};
  for (auto [is_diagonal, particle_ptr] : p.getNeighbours()) {
    auto* p2 = reinterpret_cast<Particle*>(particle_ptr);
    dvec3 rv = p2->getX() - p.getX();
    const double r = ArrayUtils::L2Norm(rv);
    dvec3 f = (k_ * (r - (is_diagonal ? sr_0_ : r_0_))) * ((1.0 / r) * rv);
    force_acc = force_acc + f;
  }
  return force_acc;
}
