//
// Created by jkr on 12/25/24.
//

#include "TruncatedLennardJones.h"

#include "utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

dvec3 TruncatedLennardJones::directionalForce(Particle& p1,
                                              Particle& p2) const {
  const dvec3 rv = p2.getX() - p1.getX();
  const double r = ArrayUtils::L2Norm(rv);
  const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
  if (r >= sigma * 1.22462048309) {
    return {0, 0, 0};
  }
  const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
  const double sr = sigma / r;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * epsilon * (sr6 - 2 * sr12) / std::pow(r, 2);

  const dvec3 force = force_magnitude * rv;

  return force;
}
