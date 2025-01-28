//
// Created by jkr on 1/28/25.
//
#include <gtest/gtest.h>

#include <vector>

#include "../src/forces/IndexForce.h"
#include "../src/utils/ArrayUtils.h"
#include "debug/debug_print.h"
#include "testUtil.h"

/*
 * TruncatedLennardJones with distance >= sigma * c
 */
TEST(IndexForce, single_index) {
  const Particle p({1, 0, 0}, {0, 0, 0}, 1, 5, 1);
  const Particle q({0, 1, 0}, {0, 0, 0}, 1, 5, 1);
  const Particle r({0, 0, 1}, {0, 0, 0}, 1, 5, 1);
  const Particle s({0, 0, 0}, {0, 0, 0}, 1, 5, 1);
  const IndexForce index_force{{1}, 10, dvec3{0.8, 0.0, 0.0}};
  double time = 5;
  std::vector<Particle> particles = {p, q, r, s};
  for (Particle &particle : particles) {
    DVEC3_NEAR(particle.getF(), {0, 0, 0}, "force not zero");
  }
  for (Particle &particle : particles) {
    for (const auto index : index_force.getIndeces()) {
      if (particle.getId() == index) {
        particle.addF(index_force.applyForce(particle, time));
      }
    }
  }
  for (Particle &particle : particles) {
    for (const auto index : index_force.getIndeces()) {
      if (particle.getId() == index) {
        DVEC3_NEAR(particle.getF(), {0.8, 0, 0}, "force");
      } else {
        DVEC3_NEAR(particle.getF(), {0, 0, 0}, "force");
      }
    }
  }
}