//
// Created by jkr on 12/23/24.
//

#include "HarmonicForce.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

dvec3 HarmonicForce::applyForce(const Particle& p) const {
  //SpdWrapper::get()->info("HarmonicForce::applyForce() on {} neighbours",
  //                         p.getNeighbours().size());
  dvec3 force_acc = {0.0, 0.0, 0.0};
  for (const auto& [fst, snd_ptr] : p.getNeighbours()) {
    if (snd_ptr) {
      const auto snd = snd_ptr;
      // SpdWrapper::get()->info(
      //     "Harmonic force between {} Particle {} and neighbour particle {}",
      //     fst ? "diagonal" : "non-diagonal", p.getId(), snd->getId());
      dvec3 rv = snd->getX() - p.getX();

      const double r = ArrayUtils::L2Norm(rv);
      // SpdWrapper::get()->info("With distance {}", r);
      dvec3 f = (k * (r - (fst ? sr_0 : r_0))) * ((1.0 / r) * rv);
      InfoVec("---", f);
      force_acc = force_acc + f;
    } else {
      SpdWrapper::get()->error("Nullptr found");
    }
  }
  InfoVec("", force_acc);
  return force_acc;
}
