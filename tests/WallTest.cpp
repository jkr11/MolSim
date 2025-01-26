//
// Created by maximilian on 22.01.25.
//
#include <gtest/gtest.h>

#include "calc/VerletIntegrator.h"
#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/LennardJones.h"
#include "forces/SingularGravity.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"

/**
 * tests that wall particles do not move
 */
TEST(Wall, immovable) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({1, 1, 1}, {1, 1, 1}, 1, 1, 1, -1);
  container.addParticle(p);
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  interactive_forces.push_back(std::make_unique<LennardJones>());

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  singular_forces.push_back(std::make_unique<SingularGravity>(1, 1));

  std::vector<std::unique_ptr<IndexForce>> index_forces;

  VerletIntegrator v(interactive_forces, singular_forces, index_forces, 0.005);
  for (int i = 0; i < 3; i++) {
    v.step(container);
  }

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "Particle Moved", 1e-8);
  });
}

TEST(Wall, excludedFromThermostat) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const Particle wall({1, 1, 1}, {1, 1, 1}, 1, 1, 1, -1);
  const Particle q({2, 2, 2}, {2, 2, 2}, 1, 1, 1, 1);

  container.addParticle(wall);
  container.addParticle(q);
  EXPECT_EQ(container.size(), 2) << "Number of Particles is not 0";
}