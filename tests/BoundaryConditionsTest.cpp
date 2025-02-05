//
// Created by maximilian on 29.11.24.
//
#include <gtest/gtest.h>

#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"

/**
 * uses an API for testing since the including part does not work for inlined
 * functions
 *
 */

/**
 * all 9 cells IN 2D should have the right declaration of type of cell
 * z-axis is not tested as it is currently not a part of the simulation
 */
TEST(BoundaryConditions, Precalculations) {
  LinkedCellsContainer container({.domain = {3, 3, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  // coordinates in the 9 different cells
  dvec3 center_center = {1, 1,
                         1};  // boundary cell and only cell of real simulation
  dvec3 center_top = {1, 4, 1};  // halo
  dvec3 center_bottom = {1, -2, 1};

  dvec3 right_top = {4, 4, 1};
  dvec3 right_center = {4, 1, 1};
  dvec3 right_bottom = {-2, -2, 1};

  dvec3 left_top = {-2, 1, 1};
  dvec3 left_center = {-2, -2, 1};
  dvec3 left_bottom = {-2, -2, 1};

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(center_center)),
            true);
  EXPECT_EQ(container.isHaloTesting(
                container.dvec3ToCellIndex_testing(center_center)),
            false);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(center_top)),
            false);
  EXPECT_EQ(
      container.isHaloTesting(container.dvec3ToCellIndex_testing(center_top)),
      true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(center_bottom)),
            false);
  EXPECT_EQ(container.isHaloTesting(
                container.dvec3ToCellIndex_testing(center_bottom)),
            true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(right_top)),
            false);
  EXPECT_EQ(
      container.isHaloTesting(container.dvec3ToCellIndex_testing(right_top)),
      true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(right_center)),
            false);
  EXPECT_EQ(container.isHaloTesting(
                container.dvec3ToCellIndex_testing(right_center)),
            true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(right_bottom)),
            false);
  EXPECT_EQ(container.isHaloTesting(
                container.dvec3ToCellIndex_testing(right_bottom)),
            true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(left_top)),
            false);
  EXPECT_EQ(
      container.isHaloTesting(container.dvec3ToCellIndex_testing(left_top)),
      true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(left_center)),
            false);
  EXPECT_EQ(
      container.isHaloTesting(container.dvec3ToCellIndex_testing(left_center)),
      true);

  EXPECT_EQ(container.isBoundaryTesting(
                container.dvec3ToCellIndex_testing(left_bottom)),
            false);
  EXPECT_EQ(
      container.isHaloTesting(container.dvec3ToCellIndex_testing(left_bottom)),
      true);
}

/**
 * if nothing can be done, nothing should be done for Outflow
 */
TEST(BoundaryConditions, Idempotence_Outflow) {
  LinkedCellsContainer container({.domain = {3, 3, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  EXPECT_EQ(container.size(), 0) << "Number of Particles is not 0";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Number of Particles is not 0";

  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";
  container.imposeInvariant();
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";
}

/**
 * test xhigh outflow
 */
TEST(BoundaryConditions, xhigh_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({4, 1, 1}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 * test xlow outflow
 */
TEST(BoundaryConditions, xlow_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({-1, 1, 1}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 * test yhigh outflow
 */
TEST(BoundaryConditions, yhigh_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({1, 4, 1}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 * test ylow outflow
 */
TEST(BoundaryConditions, ylow_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({1, -1, 1}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 * test zhigh outflow
 */
TEST(BoundaryConditions, zhigh_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({1, 1, 4}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 * test zlow outflow
 */
TEST(BoundaryConditions, zlow_Outflow) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  container.singleIterator([this](Particle& p) { p.setX({1, 1, -1}); });
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  container.imposeInvariant();
  EXPECT_EQ(container.size(), 0) << "Particle was not deleted";
}

/**
 *  test xlow Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, xlow_Reflective) {
  LinkedCellsContainer container(
      {.domain = {90, 3, 3},  // better safe than sorry
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({0.7, 1, 1}, {-1, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
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
    DVEC3_NEAR(p.getV(), {1.0, 0.0, 0.0}, "Violated the law of conservation of energy", 1e-5);
  });
}

/**
 *  test xhigh Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, xhigh_Reflective) {
  LinkedCellsContainer container(
      {.domain = {90, 3, 3},  // better safe than sorry
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({89.3, 1, 1}, {1, 0, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
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
    DVEC3_NEAR(p.getV(), {-1.0, 0.0, 0.0}, "Violated the law of conservation of energy", 1e-5);
  });
}

/**
 *  test ylow Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, ylow_Reflective) {
  LinkedCellsContainer container(
      {.domain = {3, 90, 3},  // better safe than sorry
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({1, 0.7, 1}, {0, -1, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
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
    DVEC3_NEAR(p.getV(), {0.0, 1.0, 0.0}, "Violated the law of conservation of energy", 1e-5);
  });
}

/**
 *  test yhigh Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, yhigh_Reflective) {
  LinkedCellsContainer container(
      {.domain = {3, 90, 3},  // better safe than sorry
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({1, 89.3, 1}, {0, 1, 0}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
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
    DVEC3_NEAR(p.getV(), {0.0, -1.0, 0.0}, "Violated the law of conservation of energy", 1e-5);
  });
}

/**
 *  test zlow Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, zlow_Reflective) {
  LinkedCellsContainer container(
      {.domain = {3, 90, 3},  // better safe than sorry
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
       }});

  Particle p({1, 1, 0.7}, {0, 0, -1}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
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
    DVEC3_NEAR(p.getV(), {0.0, 0, 1.0}, "Violated the law of conservation of energy", 1e-5);
  });
}

/**
 *  test zhigh Reflective boundary and that no energy is gained
 */
TEST(BoundaryConditions, zhigh_Reflective) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Reflective,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  Particle p({1, 1, 2.3}, {0, 0, 1}, 1, 1, 1);
  container.addParticles({p});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 0";

  // simulate 10.000 steps for a specific delta t to assure that it turned around
  for (int i = 0; i < 10000; i++) {
    double delta_t = 0.00005;
    container.singleIterator([this, delta_t](Particle& particle) {
      const dvec3 new_x = particle.getX() + delta_t * particle.getV() +
                          (delta_t * delta_t / (2 * particle.getM())) * (particle.getF());
      particle.setX(new_x);
    });

    container.singleIterator([](Particle& particle) { particle.updateForceInTime(); });

    container.imposeInvariant();

    container.singleIterator([this, delta_t](Particle& particle) {
      const dvec3 new_v =
          particle.getV() + (delta_t / (2 * particle.getM()) * (particle.getOldF() + particle.getF()));
      particle.setV(new_v);
    });
  }

  container.singleIterator([this](const Particle& particle) {
    DVEC3_NEAR(particle.getV(), {0.0, 0.0, -1.0}, "Violated the law of conservation of energy", 1e-5);
  });
}