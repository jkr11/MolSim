//
// Created by jkr on 12/23/24.
//

#ifndef HARMONICFORCE_H
#define HARMONICFORCE_H
#include <cmath>

#include "SingularForce.h"

class HarmonicForce : public SingularForce {
 private:
  double k{};
  double r_0{};
  double sr_0 = std::sqrt(2) * r_0;

 public:
  explicit HarmonicForce(const double k, const double r_0) : k(k), r_0(r_0) {}

  [[nodiscard]] dvec3 applyForce(const Particle& p) const override;
};

#endif  // HARMONICFORCE_H
