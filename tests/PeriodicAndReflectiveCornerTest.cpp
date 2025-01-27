//
// Created by maximilian on 22.01.25.
//
#include <gtest/gtest.h>

#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/LennardJones.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"

/**
 * tests movement through two seperate boundary conditions
 */
TEST(PeriodicAndReflective, XPeriodicYReflective) {
  LinkedCellsContainer container(
      {.domain = {9, 9, 9},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 8, 1}, {2, 2, 0}, 1, 5.0, 1);  // in bb

  container.addParticle(one);
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  // simulate 10.000 steps for a specific delta t to assure that it turned
  // around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.0001;
    container.singleIterator([this, delta_t](Particle& p) {
      const dvec3 new_x = p.getX() + delta_t * p.getV() +
                          (delta_t * delta_t / (2 * p.getM())) * (p.getF());
      p.setX(new_x);
    });

    container.singleIterator([](Particle& p) { p.updateForceInTime(); });

    container.imposeInvariant();

    container.singleIterator([this, delta_t](Particle& p) {
      const dvec3 new_v =
          p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
      p.setV(new_v);
    });
  }

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {2.0, -2.0, 0}, "wrong velocity", 1e-5);
    ASSERT_NEAR(p.getX()[0], 1, 1e-5);
    ASSERT_NEAR(p.getX()[2], 1, 1e-5);
    // y coordinate can be everything
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 2, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 2, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}