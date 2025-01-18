#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/ParticleContainer.h"
#include "../src/forces/LennardJones.h"
#include "testUtil.h"
#include "../src/utils/ArrayUtils.h"

/*
 * LennardJones directional force, arbitrary example 1
 */
TEST(LennardJones, directionalForce1) {
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 5, 1);
  Particle q({0, 1, 0}, {0, 0, 0}, 1, 5, 1);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-5.625, 5.625, 0}, "Directional force wrong", 1e-5);
}

/*
 * LennardJones directional force, arbitrary example 2
 */
TEST(LennardJones, directionalForce2) {
  Particle p({1, 2, 3}, {0, 0, 0}, 7, 5, 1);
  Particle q({4, 5, 6}, {0, 0, 0}, 10, 5, 1);
  LennardJones lj;
  dvec3 f = lj.directionalForce(p, q);
  std::cout << f << std::endl;
  DVEC3_NEAR(f, {0.000677, 0.000677, 0.000677},
             "Directional force wrong", 1e-5);
}

/*
 * LennardJones directional force, arbitrary example 3
 */
/* removed because this never happens in an actual simulation
TEST(LennardJones, directionalForce3) {
  Particle p({0, 0, 0}, {0, 0, 0}, 0, 0, 0);
  Particle q({0, 0, 0}, {0, 0, 0}, 0, 0, 0);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {0.0, 0.0, 0.0}, "Directional force wrong", 1e-5);
}
*/

/*
 * LennardJones directional force, arbitrary example 4
 */
TEST(LennardJones, directionalForce4) {
  Particle p({0.1, 1.0, 0.7}, {0, 0, 0}, 0.5, 5, 1.0);
  Particle q({0.6, 0.8, 0.3}, {0, 0, 0}, 6.4, 5, 1.0);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-30650.75270, 12260.30108, 24520.60216},
             "Directional force wrong", 1e-5);
}