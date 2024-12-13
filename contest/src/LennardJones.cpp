//
// Created by mcarn on 12/11/24.
//

#include <cmath>

#include "../../src/utils/ArrayUtils.h"
#include "LennardJones.h"

dvec3 LennardJones::directionalForce(dvec3 pos1, dvec3 pos2, double sigma1, double sigma2, double eps1, double eps2) const {
  const dvec3 differenceVector = pos2 - pos1;
  const double distance = ArrayUtils::L2Norm(differenceVector);

  const double sigma = (sigma1 + sigma2) / 2;
  const double epsilon = std::sqrt(eps1 * eps2);

  const double sr = sigma / distance;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr, 12);

  const double force_magnitude = 24 * epsilon * (sr6 - 2 * sr12) / std::pow(distance, 2);
  return force_magnitude * differenceVector;
}

double LennardJones::simpleForce(double sigma, double eps, double distance) {
  if (distance <= 0) return 0;
  const double sr = sigma / distance;
  const double sr6 = std::pow(sr, 6);
  const double sr12 = std::pow(sr6, 2);
  const double force_magnitude =
      24 * eps * (sr6 - 2 * sr12) / (std::pow(distance, 2));
  return force_magnitude;
}