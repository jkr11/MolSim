#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/DirectSumContainer.cpp"
#include "testUtil.h"

/*
 * Add particle changes internal vector length
 */
TEST(DirectSumContainer, addParticleAndSize) {
  DirectSumContainer container;
  Particle particle;

  ASSERT_EQ(container.size(), 0)
      << "ParticleContainer particle count not 0 after init.";

  container.addParticle(particle);
  ASSERT_EQ(container.size(), 1)
      << "ParticleContainer particle count not matching after addParticle.";
}

/*
 * singleIterator iterates over each distinct particle exactly once
 */
TEST(DirectSumContainer, singleIterator) {
  DirectSumContainer container;
  Particle p1;
  Particle p2;
  Particle p3;

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);

  ASSERT_EQ(container.size(), 3) << "container particle count not "
                                    "matching after adding 3 particles.";

  std::vector<Particle> vec = {};
  container.singleIterator([&vec](Particle& p) { vec.push_back(p); });

  EXPECT_VECTOR_EQ(vec, {p1, p2, p3});
}

/*
 * pair_iterator iterates over each distinct pair of particles exactly once
 */
TEST(DirectSumContainer, pair_iterator) {
  DirectSumContainer container;
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
      [&pairs](Particle& p1, Particle& p2) { pairs.emplace_back(&p1, &p2); });

  ASSERT_EQ(pairs.size(), 3);

  ASSERT_TRUE(*pairs[0].first == *container.getParticles()[0]);
  ASSERT_TRUE(*pairs[0].second == *container.getParticles()[1]);

  ASSERT_TRUE(*pairs[1].first == *container.getParticles()[0]);
  ASSERT_TRUE(*pairs[1].second == *container.getParticles()[2]);

  ASSERT_TRUE(*pairs[2].first == *container.getParticles()[1]);
  ASSERT_TRUE(*pairs[2].second == *container.getParticles()[2]);
}