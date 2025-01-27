//
// Created by jkr on 12/25/24.
//

#include "TruncatedLennardJones.h"

#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

dvec3 TruncatedLennardJones::directionalForce(Particle& p1,
                                              Particle& p2) const {
  /*for (auto [diag, ppt] : p1.getNeighbours()) {
    auto p3 = reinterpret_cast<Particle*>(ppt);
    if (p3->getId() == p2.getId()) {
      return {0, 0, 0};
    }
  }*/
  INFO("Entered LennardJonesForce")
  const dvec3 rv = p2.getX() - p1.getX();
  const double r = ArrayUtils::L2Norm(rv);
  const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
  INFO_FMT("Distance {}", r)
  INFO_FMT("Sigma * C {}", sigma * 1.1225)
  if (r >= sigma * 1.22462048309) {
    INFO("Returning 000")
    return {0, 0, 0};
  }
  INFO("Made past loop")
  const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
  const double sr = sigma / r;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * epsilon * (sr6 - 2 * sr12) / std::pow(r, 2);

  const dvec3 force = force_magnitude * rv;

  //SpdWrapper::get()->critical("Truncated lennard jones force: {}, {}, {}",
  //                            force[0], force[1], force[2]);

  return force;
}
