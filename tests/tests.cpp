//
// Created by jkr on 11/13/24.
//
#include <gtest/gtest.h>

#include <vector>

#include "../src/calc/VerletIntegrator.h"
#include "../src/defs/ParticleContainer.h"
#include "../src/forces/LennardJones.h"
#include "testUtil.h"

/*
  All tests (should) follow this convention:
  TEST(<Class>, <Function/Method>)
*/

/*
 * Add particle changes internal vector length
 */
TEST(ParticleContainer, addParticle) {
  ParticleContainer container;
  Particle particle;

  ASSERT_EQ(container.size(), 0)
      << "ParticleContainer particle count not 0 after init.";

  container.addParticle(particle);
  ASSERT_EQ(container.size(), 1)
      << "ParticleContainer particle count not matching after addParticle.";
}

/*
 * single_iterator iterates over each distinct particle exactly once
 */
TEST(ParticleContainer, single_iterator) {
  ParticleContainer container;
  Particle p1;
  Particle p2;
  Particle p3;

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);

  ASSERT_EQ(container.size(), 3) << "ParticleContainer particle count not "
                                    "matching after adding 3 particles.";

  std::vector<Particle> vec = {};
  container.single_iterator([&vec](Particle& p) { vec.push_back(p); });

  ASSERT_VECTOR_EQ(vec, {p1, p2, p3});
}

/*
 * pair_iterator iterates over each distinct pair of particles exactly once
 */
TEST(ParticleContainer, pair_iterator) {
  ParticleContainer container;
  Particle p1;
  Particle p2;
  Particle p3;

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);

  ASSERT_EQ(container.size(), 3) << "ParticleContainer particle count not "
                                    "matching after adding 3 particles.";

  std::vector<std::pair<Particle*, Particle*>> pairs;

  container.pairIterator(
      [&pairs](Particle& p1, Particle& p2) { pairs.push_back({&p1, &p2}); });

  ASSERT_EQ(pairs.size(), 3);

  ASSERT_TRUE(*pairs[0].first == container[0]);
  ASSERT_TRUE(*pairs[0].second == container[1]);

  ASSERT_TRUE(*pairs[1].first == container[0]);
  ASSERT_TRUE(*pairs[1].second == container[2]);

  ASSERT_TRUE(*pairs[2].first == container[1]);
  ASSERT_TRUE(*pairs[2].second == container[2]);
}

/*
 * LennardJones directional forces is correct
 */
TEST(LennardJones, directionalForce) {
  Particle p({1, 0, 0}, {0, 0, 0}, 1, 1, 1);
  Particle q({0, 1, 0}, {0, 0, 0}, 1, 1, 1);
  LennardJones lj;

  dvec3 f = lj.directionalForce(p, q);
  DVEC3_NEAR(f, {-1.125f, 1.125f, 0}, "Directional force wrong", 1e-5f);
}

/*
 * Positions correct after one step
 */
TEST(VerletIntegrator, step) {
  ParticleContainer container;
  Particle p({1, 0, 0}, {1, 0, 0}, 1, 1, 1);
  LennardJones lj;
  VerletIntegrator integrator(lj, 0.01f);

  p.setF({0, 1, 0});
  container.addParticle(p);
  ASSERT_EQ(container.size(), 1);

  integrator.step(container);

  // this iterator is important, since ParticleContainer is a vector
  //  therefore it necessarily copies the particle and the refrence
  //  of our p is invalid
  container.single_iterator([](Particle& p) {
    DVEC3_NEAR(p.getX(), {1.01, 0.00005, 0}, "Position wrong.", 1e-5f);
    DVEC3_NEAR(p.getOldF(), {0, 1, 0}, "Old force wrong.", 1e-5f);
    DVEC3_NEAR(p.getF(), {0, 0, 0}, "New F wrong.", 1e-5f);
    DVEC3_NEAR(p.getV(), {1, 0.005, 0}, "Velocity wrong.", 1e-5f);
  });
}
