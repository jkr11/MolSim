//
// Created by mcarn on 12/11/24.
//

#pragma once

#include "../../src/defs/types.h"
#include "InteractiveForce.h"

class LennardJones final : public InteractiveForce {
public:
  LennardJones() = default;

  dvec3 directionalForce(dvec3 pos1, dvec3 pos2, double sigma1, double sigma2, double eps1, double eps2) const;
  static double simpleForce(double sigma, double eps, double distance);
};