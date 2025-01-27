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
  for (auto [is_diagonal, particle_ptr] : p.getNeighbours()) {
    //const auto particle = particle_ptr;
    // SpdWrapper::get()->info(
    //     "Harmonic force between  Particle {} and  {} neighbour particle {}",
    ///     p.getId(), is_diagonal ? "diagonal" : "non-diagonal",
    //     particle->getId());
    // int id = particle_ptr->getId();
    // dvec3 p2_real_x = {0.0, 0.0, 0.0};
    // for (auto particle : container.getParticles()) {
    //   if (particle->getId() == id) {
    //     p2_real_x = p2_real_x + particle->getX();
    //   }
    // }

    auto* p2 = reinterpret_cast<Particle*>(particle_ptr);
    auto p2_real_x = p2->getX();
    dvec3 rv = p2->getX() - p.getX();

    const double r = ArrayUtils::L2Norm(rv);
    // SpdWrapper::get()->info("With distance {}", r);
    dvec3 f = (k * (r - (is_diagonal ? sr_0 : r_0))) * ((1.0 / r) * rv);
    // InfoVec("---", f);
    force_acc = force_acc + f;

    if (p.getId() == 874) {
      SpdWrapper::get()->info("\t874 membrane from {}: [{}, {}, {}]",
                              p2->getId(), f[0], f[1], f[2]);
      // for (auto [diag, ref] : p.getNeighbours()) {
      //   SpdWrapper::get()->info("\t\t\t{} at [{}, {}, {}]", ref.get()->getId(), particle_ptr->getX()[0], particle_ptr->getX()[1], particle_ptr->getX()[2]);
      // }
      if (p2->getId() == 823) {
        SpdWrapper::get()->info("\t\t874: [{}, {}, {}]", p.getX()[0],
                                p.getX()[1], p.getX()[2]);
        SpdWrapper::get()->info("\t\t823: [{}, {}, {}]", p2_real_x[0],
                                p2_real_x[1], p2_real_x[2]);

        SpdWrapper::get()->info("\t\tstatus {} {}", is_diagonal, r);
      }
    }

    if (p.getId() == 823 && p2->getId() == 874) {
      SpdWrapper::get()->info("\t823 membrane from {}: [{}, {}, {}]",
                              p2->getId(), f[0], f[1], f[2]);
      SpdWrapper::get()->info("\t\t823: [{}, {}, {}]", p.getX()[0], p.getX()[1],
                              p.getX()[2]);
      SpdWrapper::get()->info("\t\t874: [{}, {}, {}]", p2_real_x[0],
                              p2_real_x[1], p2_real_x[2]);
      SpdWrapper::get()->info("\t\tstatus {} {}", is_diagonal, r);
    }
  }
  // InfoVec("", force_acc);
  return force_acc;
}
