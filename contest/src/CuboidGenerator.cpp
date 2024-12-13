//
// Created by mcarn on 12/11/24.
//

#include <iostream>

#include "CuboidGenerator.h"
#include "ParticleContainer.h"
#include "../../src/utils/ArrayUtils.h"
#include "../../src/utils/MaxwellBoltzmannDistribution.h"

void CuboidGenerator::generate(ParticleContainer &particleContainer) {
  for (int i = 0; i < dimensions[0]; i++) {
    for (int j = 0; j < dimensions[1]; j++) {
      for (int k = 0; k < dimensions[2]; k++) {
        dvec3 position = {corner[0] + i * h, corner[1] + j * h, corner[2] + k * h};
        dvec3 veclocity = initialVelocity + maxwellBoltzmannDistributedVelocity(mv, twoD ? 2 : 3);

        particleContainer.addParticle(particleContainer.getNewId(), type, position, veclocity, {0, 0, 0}, {0, 0, 0});
      }
    }
  }
}