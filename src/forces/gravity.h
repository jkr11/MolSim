#pragma once
#include "../defs/Particle.h"
#include "../utils/ArrayUtils.h"
#include "force.h"

class Gravity final : public Force {
private:
  //double G = 6.674e-11; g = 1
public:
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};
