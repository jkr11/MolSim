#include <gtest/gtest.h>

#include <vector>

#include "../src/defs/containers/LinkedCellsContainer.h"
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
    LinkedCellsContainer container({10, 20, 30}, 3);

    EXPECT_EQ(container.getCellCount()[0], 5) << "X count wrong.";
    EXPECT_EQ(container.getCellCount()[1], 8) << "Y count wrong.";
    EXPECT_EQ(container.getCellCount()[2], 12) << "Z count wrong.";

    EXPECT_NEAR(container.getCellDim()[0], 3.33333, 1e-5) << "X dim wrong.";
    EXPECT_NEAR(container.getCellDim()[1], 3.33333, 1e-5) << "Y dim wrong.";
    EXPECT_NEAR(container.getCellDim()[2], 3, 1e-5) << "Z dim wrong.";
}


/*
 * .isBoundaryIndex(...) is correct for some example indices
 * .isBoundaryVec3(...) is correct for some example indices
 */
TEST(LinkedCellsContainer, isBoudaryIndex_and_isBoudaryVec3) {
  LinkedCellsContainer container({10, 20, 30}, 10); //1x2x3 => [-1, 1] x [-1, 2] x [-1, 3]

  /*EXPECT_TRUE(container.isBoundary());
  EXPECT_TRUE(container.isBoundary());
  EXPECT_TRUE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());

  EXPECT_TRUE(container.isBoundary());
  EXPECT_TRUE(container.isBoundary());
  EXPECT_TRUE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());
  EXPECT_FALSE(container.isBoundary());*/
}

/*
 * .isHaloIndex(...) is correct for some example indices
 * .isHaloVec3(...) is correct for some example indices
 */
TEST(LinkedCellsContainer, isHaloIndex_and_isHaloVec3) {
  LinkedCellsContainer container({10, 20, 30}, 10);

  /*EXPECT_TRUE(container.isHalo());
  EXPECT_TRUE(container.isHalo());
  EXPECT_TRUE(container.isHalo());
  EXPECT_FALSE(container.isHalo());
  EXPECT_FALSE(container.isHalo());
  EXPECT_FALSE(container.isHalo());

  EXPECT_TRUE(container.isHalo());
  EXPECT_TRUE(container.isHalo());
  EXPECT_TRUE(container.isHalo());
  EXPECT_FALSE(container.isHalo());
  EXPECT_FALSE(container.isHalo());
  EXPECT_FALSE(container.isHalo());*/
}


/*
 * if new container, then container.size() == 0
 * .addParticle(...) increments .size()
 * .addParticle(...) inserts into the right cell
 * .removeParticle(...) decrements .size()
 * .removeParticle(...) removes from the right cell
 */
TEST(LinkedCellsContainer, Size_addParticle_and_removeParticle) {
  LinkedCellsContainer container({10, 10, 10}, 5);
  EXPECT_EQ(container.size(), 0) << "Freshly instantiated LinkedCellsContainer is not empty.";

  Particle p = createParticle(1, 1, 1);
  container.addParticle(p);
  EXPECT_EQ(container.size(), 1) << ".addParticle() did not increase .size() by 1.";

  //TODO: cell check

  // remove particle is not coded yet
  //container.removeParticle(p);
  //EXPECT_EQ(container.size(), 0) << ".removeParticle() did not decrease .size() by 1.";

  //TODO: cell check
}

/*
 * .singleIterator() iterates over all particles
 */
TEST(LinkedCellsContainer, singleIterator) {
  LinkedCellsContainer container({10, 10, 10}, 2);

  Particle p1 = createParticle(1, 1, 1);
  Particle p2 = createParticle(5, 1, 6);
  Particle p3 = createParticle(7, 7, 8);

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
    createParticle(1, 1, 1),
    createParticle(5, 1, 6),
    createParticle(7, 7, 8),
    createParticle(4, 3, 0),
    createParticle(5, 0, 5),
    createParticle(1, 5, 2),
    createParticle(9, 6, 4),
    createParticle(2, 1, 1),
    createParticle(3, 0, 0),
    createParticle(0, 6, 1)
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

/*
 * 
 */
TEST(LinkedCellsContainer, boundaryIterator) {
  //TODO
}

/*
 *
 */
TEST(LinkedCellsContainer, haloIterator) {
  //TODO
}