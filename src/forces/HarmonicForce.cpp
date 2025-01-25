//
// Created by jkr on 12/23/24.
//

#include "HarmonicForce.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

dvec3 HarmonicForce::applyForce(const Particle& p) const {
  // SpdWrapper::get()->info(
  //"HarmonicForce::applyForce() with k : {} and r_0 : {}", k, r_0);
  dvec3 force_acc = {0.0, 0.0, 0.0};
  for (const auto& [is_diagonal, particle_ptr] : p.getNeighbours()) {
    const auto particle = particle_ptr;
    // SpdWrapper::get()->info(
    //   "Harmonic force between  Particle {} and  {} neighbour particle {}",
    //   p.getId(), is_diagonal ? "diagonal" : "non-diagonal",
    //   particle->getId());
    dvec3 rv = particle->getX() - p.getX();

    const double r = ArrayUtils::L2Norm(rv);
    // SpdWrapper::get()->info("With distance {}", r);
    dvec3 f = (k * (r - (is_diagonal ? sr_0 : r_0))) * ((1.0 / r) * rv);
    // InfoVec("---", f);
    force_acc = force_acc + f;

    if (p.getId() == 874) {
      SpdWrapper::get()->info("\t874 membrane from {}: [{}, {}, {}]", particle->getId(), f[0], f[1], f[2]);
      if (particle->getId() == 823) {
        SpdWrapper::get()->info("\t\t874: [{}, {}, {}]", p.getX()[0], p.getX()[1], p.getX()[2]);
        SpdWrapper::get()->info("\t\t823: [{}, {}, {}]", particle->getX()[0], particle->getX()[1], particle->getX()[2]);

        SpdWrapper::get()->info("\t\tstatus {} {}", is_diagonal, r);
      }
    }

    if (p.getId() == 823 && particle->getId() == 874) {
      SpdWrapper::get()->info("\t823 membrane from {}: [{}, {}, {}]", particle->getId(), f[0], f[1], f[2]);
      SpdWrapper::get()->info("\t\t823: [{}, {}, {}]", p.getX()[0], p.getX()[1], p.getX()[2]);
      SpdWrapper::get()->info("\t\t874: [{}, {}, {}]", particle->getX()[0], particle->getX()[1], particle->getX()[2]);
      SpdWrapper::get()->info("\t\tstatus {} {}", is_diagonal, r);
    }
  }
  // InfoVec("", force_acc);
  return force_acc;
}
