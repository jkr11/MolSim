//
// Created by mcarn on 12/11/24.
//

#pragma once

#include "../../src/defs/types.h"

class SingularForce {
public:
  SingularForce() = default;
  virtual ~SingularForce() = default;

  [[nodiscard]] virtual dvec3 applyForce() const = 0;
};
