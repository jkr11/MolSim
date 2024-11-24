//
// Created by jkr on 11/24/24.
//

#ifndef SINGULARGRAVITY_H
#define SINGULARGRAVITY_H
#include "SingularForce.h"

class SingularGravity final : public SingularForce {
 private:
  double g{};

 public:
  explicit SingularGravity(double g);

  // TODO: do we need to do this modular over the axis? or is it only in y
  // direction?
  [[nodiscard]] dvec3 applyForce(const Particle &p) const override;
};

#endif  // SINGULARGRAVITY_H
