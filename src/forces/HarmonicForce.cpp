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
#if 0
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
#endif
  return force_acc;
}
