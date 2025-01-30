//
// Created by jkr on 1/28/25.
//
#include <gtest/gtest.h>

#include <vector>

#include "../src/forces/IndexForce.h"
#include "../src/utils/ArrayUtils.h"
#include "debug/debug_print.h"
#include "forces/SingularGravity.h"
#include "testUtil.h"

/**
 * Tests singular gravity correctness in x direction
 */
TEST(SingularGravity, x_direction) {
  SingularGravity sg(-12.44, 0);
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 1, 1);
  p.addF(sg.applyForce(p));
  DVEC3_NEAR(p.getF(), {-12.44, 0, 0}, "Singular gravity not equal");
}

/**
 * Tests singular gravity correctness in y direction
 */
TEST(SingularGravity, y_direction) {
  SingularGravity sg(-12.44, 1);
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 1, 1);
  p.addF(sg.applyForce(p));
  DVEC3_NEAR(p.getF(), {0, -12.44, 0.0}, "Singular Gravity not equal");
}

/**
 * Tests singular gravity correctness in z direction
 */
TEST(SingularGravity, z_direction) {
  SingularGravity sg(-12.44, 2);
  Particle p({0, 0, 0}, {0, 0, 0}, 1, 1, 1);
  p.addF(sg.applyForce(p));
  DVEC3_NEAR(p.getF(), {0, 0, -12.44}, "Singular Gravity not equal");
}

/**
 * Tests singular gravity correctness for different mass
 */
TEST(SingularGravity, mass_multiplier) {
  SingularGravity sg(-12.44, 1);
  Particle p({0, 0, 0}, {0, 0, 0}, 2, 1, 1);
  p.addF(sg.applyForce(p));
  DVEC3_NEAR(p.getF(), {0, -24.88, 0.0}, "Singular Gravity not equal");
}
