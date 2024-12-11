//
// Created by mcarn on 12/11/24.
//

#pragma once

#include "../../src/defs/types.h"
#include "ParticleContainer.h"


class CuboidGenerator {
 private:
  dvec3 corner;
  ivec3 dimensions;
  double h;
  const dvec3 initialVelocity;
  double mv;
  const int type;
  const bool twoD;

 public:
  CuboidGenerator(const dvec3 &corner, const ivec3 &dimensions,
                    double h, const dvec3 &initialVelocity,
                    double mv, int type, bool twoD)
        : corner(corner),
          dimensions(dimensions),
          h(h),
          initialVelocity(initialVelocity),
          mv(mv),
          type(type),
          twoD(twoD) { }

  void generate(ParticleContainer &particleContainer);
};
