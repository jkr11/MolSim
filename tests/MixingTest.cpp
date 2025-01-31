#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/ParticleContainer.h"
#include "../src/forces/LennardJones.h"
#include "testUtil.h"
#include "../src/utils/ArrayUtils.h"
#include "../src/forces/TruncatedLennardJones.h"


/*
 * LennardJones directional force, mixing 1
 */
TEST(LennardJones, mixing) {
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 1, 2);
  Particle q({2, 2, 2}, {0, 0, 0}, 1, 3, 4);
  LennardJones lj;
  dvec3 f = lj.directionalForce(p, q);
  std::cout << f << std::endl;
  DVEC3_NEAR(f, {0.456693, 0.456693, 0.456693},
             "Directional force wrong", 1e-5);
}

/*
 * TruncatedLennardJones directional force, mixing 1
 */
TEST(TruncatedLennardJones, mixing1) {
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 1, 2);
  Particle q({2, 2, 2}, {0, 0, 0}, 1, 3, 4);
  const TruncatedLennardJones lj;
  // cutoff here is 1.1225 * sigma = 1.1225 and dist is sqrt(2) so should be 0
  const dvec3 f = lj.directionalForce(p, q);
  std::cout << f << std::endl;
  DVEC3_NEAR(f, {0.456693, 0.456693, 0.456693}, "Directional force wrong", 1e-5);
}