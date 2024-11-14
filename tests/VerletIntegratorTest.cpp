#include <gtest/gtest.h>

#include <vector>

#include "../src/calc/VerletIntegrator.h"
#include "../src/defs/ParticleContainer.h"
#include "../src/forces/LennardJones.h"
#include "testUtil.h"

/*
 * Positions correct after one step, arbitrary example 1
 */
TEST(VerletIntegrator, step1) {
  ParticleContainer container;
  Particle p({1, 0, 0}, {1, 0, 0}, 1, 1, 1);
  LennardJones lj;
  VerletIntegrator integrator(lj, 0.01f);

  p.setF({0, 1, 0});
  container.addParticle(p);
  ASSERT_EQ(container.size(), 1);

  integrator.step(container);

  p = container.getParticles()[0];
  DVEC3_NEAR(p.getX(), {1.01, 0.00005, 0}, "Position wrong.", 1e-5f);
  DVEC3_NEAR(p.getOldF(), {0, 1, 0}, "Old force wrong.", 1e-5f);
  DVEC3_NEAR(p.getF(), {0, 0, 0}, "New F wrong.", 1e-5f);
  DVEC3_NEAR(p.getV(), {1, 0.005, 0}, "Velocity wrong.", 1e-5f);
}

/*
 * Positions correct after one step, arbitrary example 2
 */
TEST(VerletIntegrator, step2) {
  ParticleContainer container;
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 1, 1);
  Particle q({0, 1, 0}, {0, 0, 0}, 1, 1, 1);
  LennardJones lj;
  VerletIntegrator integrator(lj, 0.01f);

  p.setF({0, 1, 0});
  container.addParticle(p);
  container.addParticle(q);
  ASSERT_EQ(container.size(), 2);

  integrator.step(container);

  p = container.getParticles()[0];
  DVEC3_NEAR(p.getX(), {1, 0.00005, 0}, "Position of p wrong.", 1e-5f);
  DVEC3_NEAR(p.getOldF(), {0, 1, 0}, "Old force of p wrong.", 1e-5f);
  DVEC3_NEAR(p.getF(), {-1.12516, 1.12511, 0}, "New F of p wrong.", 1e-5f);
  DVEC3_NEAR(p.getV(), {-0.005625, 0.01062, 0}, "Velocity of p wrong.", 1e-5f);

  q = container.getParticles()[1];
  DVEC3_NEAR(q.getX(), {0, 1, 0}, "Position of q wrong.", 1e-5f);
  DVEC3_NEAR(q.getOldF(), {0, 0, 0}, "Old force of q wrong.", 1e-5f);
  DVEC3_NEAR(q.getF(), {1.12516, -1.12511, 0}, "New F of q wrong.", 1e-5f);
  DVEC3_NEAR(q.getV(), {0.00562, -0.00562, 0}, "Velocity of q wrong.", 1e-5f);
}