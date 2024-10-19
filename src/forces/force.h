#include "../defs/Particle.h"
#pragma once
class Force {
public:
  virtual ~Force() = default;

  virtual dvec3 directionalForce(Particle& p1, Particle& p2) const = 0;
};
