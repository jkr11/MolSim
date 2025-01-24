//
// Created by jkr on 12/22/24.
//
#include <gtest/gtest.h>

#include <io/file/in/xml/input.cxx>
#include <io/file/in/xml/input.hxx>

#include "defs/Generators/MembraneGenerator.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/HarmonicForce.h"
#include "utils/ArrayUtils.h"

/**
 * Tests that each particle has the correct number of neighbours for the full
 * 3x3 case
 */
TEST(Membrane, 3x3x3) {
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
TEST(Membrane, LC3x3x3) {
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

  std::vector<std::unique_ptr<SingularForce>> forces;
  forces.push_back(std::make_unique<HarmonicForce>(300, 2.2));

  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0,
                                       {0, 0, 0}, 0.0, 5.0, 1.0, 1, false, {});
  membrane_generator.generate(particles);
  for (Particle& p : particles) {
    SpdWrapper::get()->info("Neighbours: {}", p.getNeighbours().size());
    for (auto& n : p.getNeighbours()) {
      if (!n.second) {
        SpdWrapper::get()->info("11111111111111 Expired pointer ---------");
      } else {
        SpdWrapper::get()->info("+++++++++++++++++++++++++++");
      }
    }
  }
  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.imposeInvariant();
  for (auto& cell : linked_cells.getCells()) {
    for (Particle& particle : cell) {
      SpdWrapper::get()->info("Number of neighbours {}",
                              particle.getNeighbours().size());
      for (auto& n : particle.getNeighbours()) {
        if (!n.second) {
          SpdWrapper::get()->info("Expired pointer ---------");
        } else {
          SpdWrapper::get()->info("+++++++++++++++++++++++++++");
        }
      }
    }
  }
  double t = 0;
  linked_cells.singleIterator([&](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (auto& force : forces) {
      SpdWrapper::get()->info("Apply force");
      f = f + force->applyForce(p);
    }
    p.addF(f);
  });
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

  std::vector<std::unique_ptr<SingularForce>> forces;
  forces.push_back(std::make_unique<HarmonicForce>(300, 2.2));

  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 1}, 1.1225, 1.0,
                                       {0, 0, 0}, 0.0, 5.0, 1.0, 1, true, {});
  membrane_generator.generate(particles);
  for (Particle& p : particles) {
    SpdWrapper::get()->info("Neighbours: {}", p.getNeighbours().size());
    for (auto& n : p.getNeighbours()) {
      if (!n.second) {
        SpdWrapper::get()->info("11111111111111 Expired pointer ---------");
      } else {
        SpdWrapper::get()->info("+++++++++++++++++++++++++++");
      }
    }
  }
  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.imposeInvariant();
  for (auto& cell : linked_cells.getCells()) {
    for (Particle& particle : cell) {
      SpdWrapper::get()->info("Number of neighbours {}",
                              particle.getNeighbours().size());
      for (auto& n : particle.getNeighbours()) {
        if (!n.second) {
          SpdWrapper::get()->info("Expired pointer ---------");
        } else {
          SpdWrapper::get()->info("+++++++++++++++++++++++++++");
        }
      }
    }
  }
  double t = 0;
  linked_cells.singleIterator([&](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (auto& force : forces) {
      SpdWrapper::get()->info("Apply force");
      f = f + force->applyForce(p);
    }
    p.addF(f);
  });
}