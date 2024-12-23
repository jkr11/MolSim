//
// Created by jkr on 12/22/24.
//
#include <gtest/gtest.h>

#include "defs/Generators/CuboidGenerator.h"

/**
 * Tests that each particle has the correct number of neighbours for the full
 * 3x3 case
 */
TEST(Membrane, 3x3) {
  std::vector<Particle> particles;
  CuboidGenerator cuboid_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0, {0, 0, 0},
                                   0.0, 5.0, 1.0, 1, false);
  cuboid_generator.generate(particles);

  ASSERT_EQ(particles[13].getNeighbours().size(), 26);
  ASSERT_EQ(particles[13].getNeighbours()[0].first, true);
  ASSERT_EQ(particles[13].getNeighbours()[1].first, true);
  ASSERT_EQ(particles[13].getNeighbours()[4].first, false);

  for (auto& particle : particles) {
    std::cout << "Address of particle: " << &particle << std::endl;
    std::cout << "Number of neighbours: " << particle.getNeighbours().size()
              << std::endl;
  }
}