//
// Created by jkr on 10/25/24.
//
#pragma once
#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "force.h"

class LennardJonesForce : public Force {
public:
  LennardJonesForce() = default;

  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};
#endif //LENNARDJONES_H
