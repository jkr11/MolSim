#include "../src/forces/TruncatedLennardJones.h"

#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/ParticleContainer.h"
#include "../src/utils/ArrayUtils.h"
#include "debug/debug_print.h"
#include "testUtil.h"

/**
 * TruncatedLennardJones with distance >= sigma * c
 */
TEST(TruncatedLennardJones, attractive_part) {
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 5, 1);
  Particle q({0, 1, 0}, {0, 0, 0}, 1, 5, 1);
  const TruncatedLennardJones lj;
  // cutoff here is 1.1225 * sigma = 1.1225 and dist is sqrt(2) so should be 0
  const dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {0, 0, 0}, "Directional force wrong", 1e-5);
}

/**
 * TruncatedLennardJones with a distance < sigma * c
 */
TEST(TruncatedLennardJones, repulsive_part) {
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 5, 1);
  Particle q({2, 0, 0}, {0, 0, 0}, 1, 5, 1);
  const TruncatedLennardJones lj;
  // cutoff here is 1.1225 * sigma = 1.1225 and dist is 1 so should be > 0
  const dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-120, 0, 0}, "Directional force wrong", 1e-5);
}

/**
 * TruncatedLennardJones at almost the cutoff
 */
TEST(TruncatedLennardJones, cutoff_part) {
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 5, 1);
  Particle q({1.1224, 0, 0}, {0, 0, 0}, 1, 5, 1);
  const TruncatedLennardJones lj;
  // cutoff here is 1.1225 * sigma = 1.1225 and dist is 1.1224 so should be > 0
  const dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-0.01774,0,0}, "Directional force wrong", 1e-5);
}