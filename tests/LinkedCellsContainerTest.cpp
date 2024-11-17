#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/LinkedCellsContainer.h"
#include "testUtil.h"

/*
 * dummy test 
 */
TEST(LinkedCellsContainer, dummy) {
    LinkedCellsContainer container({10, 10, 10}, 3);

    //EXPECT_EQ(1, 0);
}

TEST(LinkedCellsContainer, constructor) {
    LinkedCellsContainer container({10, 20, 30}, 3);

    EXPECT_EQ(container.getCellCount()[0], 5) << "X count wrong.";
    EXPECT_EQ(container.getCellCount()[1], 8) << "Y count wrong.";
    EXPECT_EQ(container.getCellCount()[2], 12) << "Z count wrong.";

    EXPECT_NEAR(container.getCellDim()[0], 3.33333, 1e-5) << "X dim wrong.";
    EXPECT_NEAR(container.getCellDim()[1], 3.33333, 1e-5) << "Y dim wrong.";
    EXPECT_NEAR(container.getCellDim()[2], 3, 1e-5) << "Z dim wrong.";
}

TEST(LinkedCellsContainer, addParticleAndSize) {
    LinkedCellsContainer container({10, 10, 10}, 5);
    EXPECT_EQ(container.size(), 0) << "Freshly instantiated LinkedCellsContainer is not empty.";

    Particle p({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
    container.addParticle(p);
    EXPECT_EQ(container.size(), 1) << ".addParticle() did not increase .size() by 1.";
}

TEST(LinkedCellsContainer, singleIterator) {
  LinkedCellsContainer container({10, 10, 10}, 2);

  Particle p1({1, 1, 1}, {0, 0, 0}, 1, 1, 1);
  Particle p2({5, 1, 6}, {0, 0, 0}, 1, 1, 1);
  Particle p3({7, 7, 8}, {0, 0, 0}, 1, 1, 1);

  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);

  EXPECT_EQ(container.size(), 3) << "container particle count not matching after adding 3 particles.";

  std::vector<Particle> vec = {};
  container.singleIterator([&vec](Particle& p) { 
    vec.push_back(p); 
  });

  EXPECT_EQ(vec.size(), 3) << "Single iterator traversed less particles than in the container.";

  EXPECT_TRUE(vec[0] == p1 || vec[1] == p1 || vec[2] == p1) << "Particle was not iterated over.";
  EXPECT_TRUE(vec[0] == p2 || vec[1] == p2 || vec[2] == p2) << "Particle was not iterated over.";
  EXPECT_TRUE(vec[0] == p3 || vec[1] == p3 || vec[2] == p3) << "Particle was not iterated over.";
}

Particle create_particle(double x, double y, double z) {
    return Particle({x, y, z}, {0, 0, 0}, 1, 1, 1);
}

/*
  Test pairIterator by running the O(n^2) algorithm and checking if 
  the count of pairs and the pairs themselves match
  Note: does not test if all pairs produces are distinct (only matters
        if total generated pair count is the same as in reference impl)
*/
TEST(LinkedCellsContainer, pairIterator) {
  const double cutoff = 5;
  LinkedCellsContainer container({10, 10, 10}, cutoff);

  std::array<Particle, 10> particles = {
    create_particle(1, 1, 1),
    create_particle(5, 1, 6),
    create_particle(7, 7, 8),
    create_particle(4, 3, 0),
    create_particle(5, 0, 5),
    create_particle(1, 5, 2),
    create_particle(9, 6, 4),
    create_particle(2, 1, 1),
    create_particle(3, 0, 0),
    create_particle(0, 6, 1)
  };

  for (int i = 0; i < particles.size(); i++) {
    container.addParticle(particles[i]);
  }

  EXPECT_EQ(container.size(), particles.size()) << "container particle count not matching after adding 4 particles.";

  //compute pairs using the slow method
  std::vector<std::array<Particle*, 2>> pairs = {};
  for (int i = 0; i < particles.size(); i++) {
    for (int j = i + 1; j < particles.size(); j++) {
      auto posp = (particles[i]).getX();
      auto posq = (particles[j]).getX();
      dvec3 d = {posp[0] - posq[0], posp[1] - posq[1], posp[2] - posq[2]};

      if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > cutoff * cutoff)
        continue;

      pairs.push_back({&particles[i], &particles[j]});
    }
  }

  int count = 0;
  container.pairIterator([&pairs, &count](Particle& p, Particle &q) {
    count++;

    for (auto it = pairs.begin(); it != pairs.end(); ++it) {
      if ((*(*it)[0] == p && *(*it)[1] == q) || (*(*it)[0] == q && *(*it)[1] == p))
        return;
    }

    EXPECT_TRUE(false) << "Pair Iterator produced a invalid pair";
  });

  EXPECT_EQ(count, pairs.size()) << "Pair count does not match reference implementation";
}