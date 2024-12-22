//
// Created by jkr on 12/22/24.
//
#include <gtest/gtest.h>

#include "defs/Generators/CuboidGenerator.h"

TEST(Membrane, 3x3) {
  std::vector<Particle> particles;
  CuboidGenerator cuboid_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0, {0, 0, 0},
                                   0.0, 5.0, 1.0, 1, false);
  cuboid_generator.generate(particles);

  ASSERT_EQ(particles[13].neighbours.size(), 26);

  for (auto& particle : particles) {
    std::cout << "Address of particle: " << &particle << std::endl;
    std::cout << "Number of neighbours: " << particle.neighbours.size()
              << std::endl;
  }
}