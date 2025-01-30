//
// Created by jkr on 12/25/24.
//

#ifndef TRUNCATEDLENNARDJONES_H
#define TRUNCATEDLENNARDJONES_H
#include "InteractiveForce.h"

class TruncatedLennardJones final : public InteractiveForce {
 public:
  explicit TruncatedLennardJones() = default;

  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};

#endif  // TRUNCATEDLENNARDJONES_H
