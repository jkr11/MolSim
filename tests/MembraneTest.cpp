//
// Created by jkr on 12/22/24.
//
#include <gtest/gtest.h>

#include <io/file/in/xml/input.cxx>
#include <io/file/in/xml/input.hxx>

#include "debug/debug_print.h"
#include "defs/Generators/MembraneGenerator.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/HarmonicForce.h"
#include "forces/IndexForce.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularGravity.h"
#include "forces/TruncatedLennardJones.h"
#include "testUtil.h"
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
  forces.push_back(std::make_unique<HarmonicForce>(3, 2.2));

  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 3}, 1.1225, 1.0,
                                       {0, 0, 0}, 0.0, 5.0, 1.0, 1, false, {});
  membrane_generator.generate(particles);

  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.imposeInvariant();
  for (auto& cell : linked_cells.getCells()) {
    for (auto particle : cell) {
      SpdWrapper::get()->info("Number of neighbours {}",
                              particle->getNeighbours().size());
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

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  singular_forces.push_back(std::make_unique<SingularGravity>(-0.001, 2));
  singular_forces.push_back(std::make_unique<HarmonicForce>(300, 2.2));

  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  interactive_forces.push_back(std::make_unique<TruncatedLennardJones>());

  std::vector<std::unique_ptr<IndexForce>> index_forces;
  IndexForceConfig index_force_config = {
      .ids = {5},
      .end_time = 15,
      .force_values = {0, 0, 0.8},
  };
  index_forces.push_back(std::make_unique<IndexForce>());

  MembraneGenerator membrane_generator({0, 0, 0}, {3, 3, 1}, 2.2, 1.0,
                                       {0, 0, 0}, 0.0, 1.0, 1.0, 1, true, {});
  membrane_generator.generate(particles);
  for (Particle& p : particles) {
    SpdWrapper::get()->info("Neighbours: {}", p.getNeighbours().size());
    for (const auto& [fst, snd] : p.getNeighbours()) {
      // check for expired pointer (legacy from weak_ptr impl, maybe we go back
      ASSERT_TRUE(snd);
    }
  }
  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.imposeInvariant();
  for (auto& cell : linked_cells.getCells()) {
    for (auto particle : cell) {
      SpdWrapper::get()->info("Number of neighbours {}",
                              particle->getNeighbours().size());
      for (const auto& [fst, snd] : particle->getNeighbours()) {
        ASSERT_TRUE(snd);
      }
    }
  }
  double t = 0;
  for (int i = 0; i < 5; i++) {
    linked_cells.singleIterator([&](Particle& p) {
      dvec3 f = {0, 0, 0};
      for (const auto& force : singular_forces) {
        SpdWrapper::get()->info("Apply force");
        f = f + force->applyForce(p);
      }
      InfoVec("total singular forces {}", f);
      p.addF(f);
    });
  }
}

TEST(Membrane, LC4x4) {
  constexpr LinkedCellsConfig linked_cells_config = {
      .domain = {20, 20, 20},
      .cutoff_radius = 4.0,
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

  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  // singular_forces.push_back(std::make_unique<SingularGravity>(-0.001, 2));
  singular_forces.push_back(std::make_unique<HarmonicForce>(300, 2.2));

  // std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  // interactive_forces.push_back(std::make_unique<TruncatedLennardJones>());

  // std::vector<std::unique_ptr<IndexForce>> index_forces;
  // IndexForceConfig index_force_config = {
  //    .ids = {5},
  //    .end_time = 15,
  //    .force_values = {0, 0, 0.8},
  // };
  // index_forces.push_back(std::make_unique<IndexForce>());
  std::cout << "global id counter " << std::endl;
  MembraneGenerator membrane_generator({0, 0, 0}, {4, 4, 1}, 2.5, 1.0,
                                       {0, 0, 0}, 0.0, 1.0, 1.0, 1, true, {});

  membrane_generator.generate(particles);

  ASSERT_EQ(particles[3].getId(), 3) << "Particle 3 has wrong id";
  ASSERT_EQ(particles[3].getNeighbours().size(), 3)
      << "Particle 3 has wrong amount of neighbours";

  for (int i = 0; i < 3; ++i) {
    Particle* neighbourParticle =
        reinterpret_cast<Particle*>(particles[3].getNeighbours()[i].second);
    int idN = neighbourParticle->getId();
    bool diag = particles[3].getNeighbours()[i].first;
    switch (idN) {
      case 2:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 2 has wrong directionality";
        break;
      case 6:
        EXPECT_TRUE(diag)
            << "Particle 3 neighbour with id 6 has wrong directionality";
        break;
      case 7:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      default:
        FAIL();
        break;
    }
  }

  for (int i = 0; i < 3; i++) {
    Particle* neighbourParticle =
        reinterpret_cast<Particle*>(particles[5].getNeighbours()[i].second);
    int idN = neighbourParticle->getId();

    bool diag = particles[5].getNeighbours()[i].first;
    switch (idN) {
      case 0:
        EXPECT_TRUE(diag)
            << "Particle 3 neighbour with id 2 has wrong directionality";
        break;
      case 1:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 6 has wrong directionality";
        break;
      case 2:
        EXPECT_TRUE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      case 6:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      case 10:
        EXPECT_TRUE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      case 9:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      case 8:
        EXPECT_TRUE(diag)
            << "Particle 3 neighbour with id 7 has wrong directionality";
        break;
      case 4:
        EXPECT_FALSE(diag)
            << "Particle 3 neighbour with id 4 has wrong directionality";
        break;
      default:
        FAIL();
        break;
    }
  }

  DVEC3_NEAR(particles[3].getX(), {0, 7.5, 0}, "Position is wrong", 0.00001);
  DVEC3_NEAR(particles[5].getX(), {2.5, 2.5, 0}, "Position o 5 is wrong",
             0.00001);

  dvec3 p1harmonic_force = singular_forces[0]->applyForce(particles[0]);
  // [0] -> [1] == (0,90,0)
  // [0] -> [4] == (90,0,0)
  // [0] -> [5] == (90,90,0)
  DVEC3_NEAR(p1harmonic_force, {180, 180, 0}, "Force", 1e-8);

  LinkedCellsContainer linked_cells(linked_cells_config);
  linked_cells.addParticles(particles);
  linked_cells.singleIterator(
      [&](Particle& p) { p.addF(singular_forces[0]->applyForce(p)); });
}