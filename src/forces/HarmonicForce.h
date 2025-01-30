//
// Created by jkr on 12/23/24.
//
#ifndef HARMONICFORCE_H
#define HARMONICFORCE_H
#include <cmath>

#include "SingularForce.h"

class HarmonicForce final : public SingularForce {
 private:
  int k_{};
  double r_0_{};
  double sr_0_{};

 public:
  explicit HarmonicForce(const int k, const double r_0) : k_(k), r_0_(r_0) {
    sr_0_ = std::sqrt(2) * r_0;
  }

  [[nodiscard]] dvec3 applyForce(const Particle& p) const override;
};

#endif  // HARMONICFORCE_H
