//
// Created by jkr on 10/18/24.
//
#pragma once
#include "../defs/Particle.h"
#include "../utils/ArrayUtils.h"
#include "force.h"

class Gravity final : public Force {
public:
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};
