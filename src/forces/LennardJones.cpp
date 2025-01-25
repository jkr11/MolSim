//
// Created by jkr on 10/25/24.
//
#include "LennardJones.h"

#include "utils/ArrayUtils.h"  // do not remove this even if clion copes about it
#include "utils/SpdWrapper.h"

dvec3 LennardJones::directionalForce(Particle& p1, Particle& p2) const {
  //SpdWrapper::get()->error("huhu");
  const dvec3 r = p2.getX() - p1.getX();
  const double rsquared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  const double sigmasquareddivrsquared =
      std::pow(p1.getSigma() + p2.getSigma(), 2) / (4 * rsquared);
  const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
  const double sr3 = std::pow(sigmasquareddivrsquared, 3);
  const double sr6 = 2 * std::pow(sr3, 2);
  const double force_magnitude = 24 * epsilon * (sr3 - sr6) / rsquared;

  return force_magnitude * r;
}

double LennardJones::simpleForce(const Particle& p, const double distance) {
  if (distance <= 0) return 0;
  const double sr = p.getSigma() / distance;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * p.getEpsilon() * (sr6 - 2 * sr12) / (std::pow(distance, 2));
  return force_magnitude;
}

dvec3 LennardJones::directionalForceWithOffset(const Particle& p1,
                                               const Particle& p2,
                                               const dvec3& distance) {
  const double r = ArrayUtils::L2Norm(distance);
  if (r == 0) return {0, 0, 0};
  const double sigma = (p1.getSigma() + p2.getSigma()) / 2;
  const double epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
  const double sr = sigma / r;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * epsilon * (sr6 - 2 * sr12) / std::pow(r, 2);
  return force_magnitude * distance;
}
