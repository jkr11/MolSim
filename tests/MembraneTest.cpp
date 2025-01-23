//
// Created by jkr on 12/22/24.
//
#include <gtest/gtest.h>

#include <io/file/in/xml/input.hxx>

#include "defs/Generators/MembraneGenerator.h"
#include "defs/containers/LinkedCellsContainer.h"

/**
 * Tests that each particle has the correct number of neighbours for the full
 * 3x3 case
 */
TEST(Membrane, 3x3) {
  std::vector<Particle> particles;
  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0,
                                       {0, 0, 0}, 0.0, 5.0, 1.0, 1, false, {});
  membrane_generator.generate(particles);

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

/**
 * Test that neighbours get correctly translated to the linkedcellscontainer
 * (Technically this is fixed by default constructor)
 */
TEST(Membrane, LC3x3) {
  constexpr LinkedCellsConfig linked_cells_config = {
      .domain = {10, 10, 10},
      .cutoff_radius = 3.0,
      .boundary_config =
          {
              .x_high = LinkedCellsConfig::Outflow,
              .x_low = LinkedCellsConfig::Outflow,
              .y_high = LinkedCellsConfig::Outflow,
              .y_low = LinkedCellsConfig::Outflow,
              .z_high = LinkedCellsConfig::Outflow,
              .z_low = LinkedCellsConfig::Outflow,
          },
  };
  std::vector<Particle> particles;

  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0,
                                       {0, 0, 0}, 0.0, 5.0, 1.0, 1, false, {});
  membrane_generator.generate(particles);
  for (Particle& p : particles) {
    SpdWrapper::get()->info("Neighbours: {}", p.getNeighbours().size());
  }
  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.imposeInvariant();
  for (auto& cell : linked_cells.getCells()) {
    for (Particle& particle : cell) {
      SpdWrapper::get()->info("Number of neighbours {}",
                              particle.getNeighbours().size());
    }
  }
}