//
// Created by jkr on 12/25/24.
//

#ifndef LENNARDJONESCUTOFF_H
#define LENNARDJONESCUTOFF_H
#include "InteractiveForce.h"

class LennardJonesCutoff final : public InteractiveForce {
 private:
  double r_cutoff{};

 public:
  explicit LennardJonesCutoff(const double r_cutoff) : r_cutoff(r_cutoff) {}

  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};

#endif  // LENNARDJONESCUTOFF_H
