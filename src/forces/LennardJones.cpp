//
// Created by jkr on 10/25/24.
//
#include "LennardJones.h"

#include "utils/ArrayUtils.h"  // do not remove this even if clion copes about it
#include "utils/SpdWrapper.h"

dvec3 LennardJones::directionalForce(Particle& p1, Particle& p2) const {
  const dvec3 rv = p2.getX() - p1.getX();
  const double r = ArrayUtils::L2Norm(rv);
  if (r == 0) return {0, 0, 0};
  constexpr double sigma = 1.0;
  constexpr double epsilon = 5.0;
  const double sr = sigma / r;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * epsilon * (sr6 - 2 * sr12) / std::pow(r, 2);
  return force_magnitude * rv;
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

dvec3 LennardJones::directionalForceWithOffset(Particle& p1, Particle& p2,
                                               const dvec3& distance) {
  const double r = ArrayUtils::L2Norm(distance);

  SpdWrapper::get()->info("[{}, {}, {}] -> {}", distance[0], distance[1], distance[2], r);

  if (r == 0) return {0, 0, 0};
  constexpr double sigma = 1.0;
  constexpr double epsilon = 5.0;
  const double sr = sigma / r;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * epsilon * (sr6 - 2 * sr12) / std::pow(r, 2);
  return force_magnitude * distance;
}

