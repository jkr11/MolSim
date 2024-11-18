#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/ParticleContainer.h"
#include "../src/forces/LennardJones.h"
#include "testUtil.h"

/*
 * LennardJones directional force, arbitrary example 1
 */
TEST(LennardJones, directionalForce1) {
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 1, 1);
  Particle q({0, 1, 0}, {0, 0, 0}, 1, 1, 1);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-1.125, 1.125, 0}, "Directional force wrong", 1e-5);
}

/*
 * LennardJones directional force, arbitrary example 2
 */
TEST(LennardJones, directionalForce2) {
  Particle p({1, 2, 3}, {0, 0, 0}, 7, 8, 9);
  Particle q({4, 5, 6}, {0, 0, 0}, 10, 11, 12);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-230211.96555, -230211.96555, -230211.96555},
             "Directional force wrong", 1e-5);
}

/*
 * LennardJones directional force, arbitrary example 3
 */
TEST(LennardJones, directionalForce3) {
  Particle p({0, 0, 0}, {0, 0, 0}, 0, 0, 0);
  Particle q({0, 0, 0}, {0, 0, 0}, 0, 0, 0);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {0.0, 0.0, 0.0}, "Directional force wrong", 1e-5);
}

/*
 * LennardJones directional force, arbitrary example 4
 */
TEST(LennardJones, directionalForce4) {
  Particle p({0.1, 1.0, 0.7}, {0, 0, 0}, 0.5, 0.9, 1.2);
  Particle q({0.6, 0.8, 0.3}, {0, 0, 0}, 6.4, 2.0, 1.3);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-123897.24283, 49558.89713, 99117.79426},
             "Directional force wrong", 1e-5);
}