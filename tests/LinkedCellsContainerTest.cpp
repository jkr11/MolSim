#include <gtest/gtest.h>

#include <vector>

#include "defs/Simulation.h"
#include "defs/containers/LinkedCellsContainer.cpp"
#include "defs/containers/LinkedCellsContainer.h"
#include "defs/containers/ParticleContainer.h"
#include "defs/types.h"
#include "testUtil.h"

/*
 * @brief small helper function to create particles
 * @param x x coordinate of particle
 * @param y y coordinate of particle
 * @param z z coordinate of particle
 * @return Particle
 */
Particle createParticle(double x, double y, double z) {
  return Particle({x, y, z}, {0, 0, 0}, 1, 1, 1);
}

/*
 * Container creates the right amount of cells and correct cell dimensions
 */
TEST(LinkedCellsContainer, constructor) {
  LinkedCellsConfig config = {.domain = {10, 20, 30},
                              .cutoff_radius = 3,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  EXPECT_EQ(container.getCellCount()[0], 5) << "X count wrong.";
  EXPECT_EQ(container.getCellCount()[1], 8) << "Y count wrong.";
  EXPECT_EQ(container.getCellCount()[2], 12) << "Z count wrong.";

  EXPECT_NEAR(container.getCellDim()[0], 3.33333, 1e-5) << "X dim wrong.";
  EXPECT_NEAR(container.getCellDim()[1], 3.33333, 1e-5) << "Y dim wrong.";
  EXPECT_NEAR(container.getCellDim()[2], 3, 1e-5) << "Z dim wrong.";
}

/*
 * .isBoundary(...) is correct for some examples
 */
TEST(LinkedCellsContainer, isBoudary) {
  LinkedCellsConfig config = {.domain = {20, 20, 30},
                              .cutoff_radius = 10,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(
      config);  // 2x2x3 => [-1, 2] x [-1, 2] x [-1, 3]
  EXPECT_TRUE(container.isBoundary({0, 0, 0}));
  EXPECT_TRUE(container.isBoundary({0, 1, 0}));
  EXPECT_TRUE(container.isBoundary({1, 1, 2}));
  EXPECT_FALSE(container.isBoundary({-1, 0, 0}));
  EXPECT_FALSE(container.isBoundary({0, 2, 0}));
  EXPECT_FALSE(container.isBoundary({1, 0, 3}));

  EXPECT_TRUE(container.isBoundary(26));
  EXPECT_TRUE(container.isBoundary(31));
  EXPECT_TRUE(container.isBoundary(52));
  EXPECT_FALSE(container.isBoundary(6));
  EXPECT_FALSE(container.isBoundary(36));
  EXPECT_FALSE(container.isBoundary(49));
}

/*
 * .isHalo(...) is correct for some examples
 */
TEST(LinkedCellsContainer, isHalo) {
  LinkedCellsConfig config = {.domain = {20, 20, 30},
                              .cutoff_radius = 10,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(
      config);  // 2x2x3 => [-1, 2] x [-1, 2] x [-1, 3]

  EXPECT_TRUE(container.isHalo({-1, -1, -1}));
  EXPECT_TRUE(container.isHalo({0, 2, 0}));
  EXPECT_TRUE(container.isHalo({-1, 2, 3}));
  EXPECT_FALSE(container.isHalo({0, 0, 0}));
  EXPECT_FALSE(container.isHalo({1, 1, 2}));
  EXPECT_FALSE(container.isHalo({0, 1, 0}));

  EXPECT_TRUE(container.isHalo(0));
  EXPECT_TRUE(container.isHalo(36));
  EXPECT_TRUE(container.isHalo(39));
  EXPECT_FALSE(container.isHalo(26));
  EXPECT_FALSE(container.isHalo(53));
  EXPECT_FALSE(container.isHalo(31));
}

/*
 * .cellIndexToCoord(...) works for some examples
 * .cellCoordToIndex(...) works for some exmaples
 * .isValidCellCoordinate(...) works for some examples
 */
TEST(LinkedCellsContainer,
     cellIndexToCoord_cellCoordToIndex_isValidCellCoordinate) {
  LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 3,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  EXPECT_EQ(container.cellIndexToCoord(0)[0], -1)
      << ".cellIndexToCoord(...) x coordinate wrong";
  EXPECT_EQ(container.cellIndexToCoord(0)[1], -1)
      << ".cellIndexToCoord(...) y coordinate wrong";
  EXPECT_EQ(container.cellIndexToCoord(0)[2], -1)
      << ".cellIndexToCoord(...) z coordinate wrong";

  EXPECT_EQ(container.cellIndexToCoord(69)[0], 1)
      << ".cellIndexToCoord(...) x coordinate wrong";
  EXPECT_EQ(container.cellIndexToCoord(69)[1], 2)
      << ".cellIndexToCoord(...) y coordinate wrong";
  EXPECT_EQ(container.cellIndexToCoord(69)[2], 3)
      << ".cellIndexToCoord(...) z coordinate wrong";

  EXPECT_EQ(container.cellCoordToIndex({-1, -1, -1}), 0)
      << ".cellCoordToIndex(...) index wrong";
  EXPECT_EQ(container.cellCoordToIndex({1, 2, 3}), 69)
      << ".cellCoordToIndex(...) index wrong";

  EXPECT_TRUE(container.isValidCellCoordinate({-1, -1, -1}))
      << ".isValidCellCoordinate(...) produced wrong result";
  EXPECT_TRUE(container.isValidCellCoordinate({0, 0, 1}))
      << ".isValidCellCoordinate(...) produced wrong result";

  EXPECT_FALSE(container.isValidCellCoordinate({-2, -1, -1}))
      << ".isValidCellCoordinate(...) produced wrong result";
  EXPECT_FALSE(container.isValidCellCoordinate({-1, 5, -1}))
      << ".isValidCellCoordinate(...) produced wrong result";
}

/*
 * if new container, then container.size() == 0
 * .addParticle(...) increments .size()
 * .removeParticle(...) decrements .size()
 */
TEST(LinkedCellsContainer, Size_addParticle_and_removeParticle) {
  LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 5,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);
  EXPECT_EQ(container.size(), 0)
      << "Freshly instantiated LinkedCellsContainer is not empty.";

  Particle p = createParticle(1, 1, 1);
  container.addParticle(p);
  EXPECT_EQ(container.size(), 1)
      << ".addParticle() did not increase .size() by 1.";

  Particle pr{{}};
  container.singleIterator([&pr](const Particle& q) { pr = q; });

  container.removeParticle(pr);
  EXPECT_EQ(container.size(), 0)
      << ".removeParticle() did not decrease .size() by 1.";
}

/*
 * .singleIterator() iterates over all particles
 */
TEST(LinkedCellsContainer, singleIterator) {
  LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 2,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  Particle p1 = createParticle(1, 1, 1);
  Particle p2 = createParticle(5, 1, 6);
  Particle p3 = createParticle(7, 7, 8);

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);

  EXPECT_EQ(container.size(), 3)
      << "container particle count not matching after adding 3 particles.";

  std::vector<Particle> vec = {};
  container.singleIterator([&vec](Particle& p) { vec.push_back(p); });

  EXPECT_EQ(vec.size(), 3)
      << "Single iterator traversed less particles than in the container.";

  EXPECT_TRUE(vec[0] == p1 || vec[1] == p1 || vec[2] == p1)
      << "Particle was not iterated over.";
  EXPECT_TRUE(vec[0] == p2 || vec[1] == p2 || vec[2] == p2)
      << "Particle was not iterated over.";
  EXPECT_TRUE(vec[0] == p3 || vec[1] == p3 || vec[2] == p3)
      << "Particle was not iterated over.";
}

/*
  Test pairIterator by running the O(n^2) algorithm and checking if
  the count of pairs and the pairs themselves match
  Note: does not test if all pairs produces are distinct (only matters
        if total generated pair count is the same as in reference impl)
*/
TEST(LinkedCellsContainer, pairIterator) {
  constexpr LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 5,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  std::array<Particle, 10> particles = {
      createParticle(1, 1, 1),
      createParticle(5, 1, 6),
      createParticle(7, 7, 8),
      createParticle(4, 3, 0),
      createParticle(5, 0, 5),
      createParticle(1, 5, 2),
      createParticle(9, 6, 4),
      createParticle(2, 1, 1),
      createParticle(3, 0, 0),
      createParticle(0, 6, 1)};

  for (std::size_t i = 0; i < particles.size(); i++) {
    container.addParticle(particles[i]);
  }

  EXPECT_EQ(container.size(), particles.size())
      << "container particle count not matching after adding 4 particles.";

  // compute pairs using the slow method
  std::vector<std::array<Particle*, 2>> pairs = {};
  for (std::size_t i = 0; i < particles.size(); i++) {
    for (std::size_t j = i + 1; j < particles.size(); j++) {
      auto posp = (particles[i]).getX();
      auto posq = (particles[j]).getX();
      dvec3 d = {posp[0] - posq[0], posp[1] - posq[1], posp[2] - posq[2]};

      if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] >
          config.cutoff_radius * config.cutoff_radius)
        continue;

      pairs.push_back({&particles[i], &particles[j]});
    }
  }


  std::cout << "Pairs are: " << std::endl;
  for (auto pair : pairs) {
    std::cout << pair[0]->getId() << " and " << pair[1]->getId() << std::endl;
  }
  std::cout << "-------------" << std::endl;


  int count = 0;
  container.pairIterator([&pairs, &count](const Particle& p,
                                          const Particle& q) {
    count++;
    std::cout<< "count is now " << count << std::endl;

    std::cout << p.getId() << " and " << q.getId() << std::endl;

    for (auto & pair : pairs) {
      if ((pair[0]->getId() == p.getId() && pair[1]->getId() == q.getId()) ||
          (pair[0]->getId() == q.getId() && pair[1]->getId() == p.getId()))
        return;
    }

    FAIL() << "Pair Iterator produced an invalid pair " << p.getId() << " and " << q.getId();
  });

  EXPECT_EQ(count, pairs.size())
      << "Pair count does not match reference implementation";
}

/*
 * boundaryIterator(...) goes over all particles in boundary
 */
TEST(LinkedCellsContainer, boundaryIterator) {
  LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 3,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  Particle p1 = createParticle(0, 0, 0);
  Particle p2 = createParticle(9, 0, 0);
  Particle p3 = createParticle(0, 9, 3);
  Particle p4 = createParticle(4, 4, 4);

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);
  container.addParticle(p4);

  container.boundaryIterator([&p1, &p2, &p3, &p4](Particle& p) {
    EXPECT_TRUE(p == p1 || p == p2 || p == p3 || !(p == p4));
  });
}

/*
 * haloIterator(...) goes over all particles in halo
 */
TEST(LinkedCellsContainer, haloIterator) {
  LinkedCellsConfig config = {.domain = {10, 10, 10},
                              .cutoff_radius = 5,
                              .boundary_config = {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }};

  LinkedCellsContainer container(config);

  Particle p1 = createParticle(-1, -1, -1);
  Particle p2 = createParticle(-1, 0, 0);
  Particle p3 = createParticle(0, 11, 0);
  Particle p4 = createParticle(0, 0, 0);

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);
  container.addParticle(p4);

  container.haloIterator([&p1, &p2, &p3, &p4](Particle& p) {
    EXPECT_TRUE(p == p1 || p == p2 || p == p3 || !(p == p4));
  });
}
