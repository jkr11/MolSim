//
// Created by maximilian on 10.12.24.
//
#include <gtest/gtest.h>

#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/LennardJones.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"

/**global variables for use in all tests
 *    +--+--+--+--+--+
 *  3 |ea|eb|ec|ed|ee|
 *    +--+--+--+--+--+
 *  2 |da|db|dc|dd|de|
 *    +--+--+--+--+--+
 *  1 |ca|cb|cc|cd|ce|
 *    +--+--+--+--+--+
 *  0 |ba|bb|bc|bd|be|
 *    +--+--+--+--+--+
 * -1 |aa|ab|ac|ad|ae|
 *    +--+--+--+--+--+
 *     -1  0  1  2  3
 */
ivec3 ea = {-1, 3, 0};
ivec3 eb = {0, 3, 0};
ivec3 ec = {1, 3, 0};
ivec3 ed = {2, 3, 0};
ivec3 ee = {3, 3, 0};
ivec3 da = {-1, 2, 0};
ivec3 db = {0, 2, 0};
ivec3 dc = {1, 2, 0};
ivec3 dd = {2, 2, 0};
ivec3 de = {3, 2, 0};
ivec3 ca = {-1, 2, 0};
ivec3 cb = {0, 2, 0};
ivec3 cc = {1, 2, 0};
ivec3 cd = {2, 2, 0};
ivec3 ce = {3, 2, 0};
ivec3 ba = {-1, 0, 0};
ivec3 bb = {0, 0, 0};
ivec3 bc = {1, 0, 0};
ivec3 bd = {2, 0, 0};
ivec3 be = {3, 0, 0};
ivec3 aa = {-1, -1, 0};
ivec3 ab = {0, -1, 0};
ivec3 ac = {1, -1, 0};
ivec3 ad = {2, -1, 0};
ivec3 ae = {3, -1, 0};

//[[======================== warping tests ===================================]]

/**
 * tests that warping is correct for x cells in 2D for a 3x3 cell grid
 * with both axis Periodic
 */
TEST(PeriodicBoundary, warpingBothPeriodicForX) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ee, xhigh);
  constexpr dvec3 expected_offset1 = {9, 9, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset1, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(de, xhigh);
  constexpr dvec3 expected_offset2 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ce, xhigh);
  constexpr dvec3 expected_offset3 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, cb);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(be, xhigh);
  constexpr dvec3 expected_offset4 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ae, xhigh);
  constexpr dvec3 expected_offset5 = {9, -9, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset5, "Offset wrong", 1e-5);
}

/**
 * tests that warping is correct for y cells in 2D for a 3x3 cell grid
 * with both axis Periodic
 */
TEST(PeriodicBoundary, warpingBothPeriodicForY) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ea, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency1";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(eb, yhigh);
  constexpr dvec3 expected_offset2 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency2";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ec, yhigh);
  constexpr dvec3 expected_offset3 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency3";
  EXPECT_IVEC3_EQ(new_coordinates, bc);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ed, yhigh);
  constexpr dvec3 expected_offset4 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency4";
  EXPECT_IVEC3_EQ(new_coordinates, bd);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ee, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency5";
}

/**
 * tests that warping is correct for x cells in 2D for a 3x3 cell grid
 * with x-axis Periodic
 */
TEST(PeriodicBoundary, warpingXPeriodicForX) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ee, xhigh);
  EXPECT_EQ(is_adjacent, false);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(de, xhigh);
  constexpr dvec3 expected_offset2 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ce, xhigh);
  constexpr dvec3 expected_offset3 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, cb);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(be, xhigh);
  constexpr dvec3 expected_offset4 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ae, xhigh);
  EXPECT_EQ(is_adjacent, false);
}

/**
 * tests that warping is correct for y cells in 2D for a 3x3 cell grid
 * with y-axis Periodic
 */
TEST(PeriodicBoundary, warpingYPeriodicForY) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ea, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency1";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(eb, yhigh);
  constexpr dvec3 expected_offset2 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency2";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ec, yhigh);
  constexpr dvec3 expected_offset3 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency3";
  EXPECT_IVEC3_EQ(new_coordinates, bc);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ed, yhigh);
  constexpr dvec3 expected_offset4 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency4";
  EXPECT_IVEC3_EQ(new_coordinates, bd);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflectiveWarpAroundTesting(ee, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency5";
}

//[[======================== force tests =====================================]]

/**
 * tests that the correct force ist calculated for a periodic boundary in the x
 * dimension
 */
TEST(PeriodicBoundaryForce, offsetX) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({8, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({-1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting(be, xhigh);
  constexpr dvec3 expected_offset = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  ASSERT_EQ(LennardJones::directionalForceWithOffset(
                *p, *q, accounted_particle_distance),
            f.directionalForce(*test, *q));
}

/**
 * tests that the correct force ist calculated for a periodic boundary in the y
 * dimension
 */
TEST(PeriodicBoundaryForce, offsetY) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({1, 8, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({1, -1, 1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting(eb, yhigh);
  constexpr dvec3 expected_offset = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  ASSERT_EQ(LennardJones::directionalForceWithOffset(
                *p, *q, accounted_particle_distance),
            f.directionalForce(*test, *q));
}

/**
 * tests that the correct force ist calculated for a periodic boundary in the z
 * dimension
 */
TEST(PeriodicBoundaryForce, offsetZ) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 9},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({1, 1, 8}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({1, 1, -1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting({0, 0, 3}, zhigh);
  constexpr dvec3 expected_offset = {0, 0, 9};

  SpdWrapper::get()->info("{}, [{}, {}, {}], [{}, {}, {}]", is_adjacent,
                          new_coordinates[0], new_coordinates[1],
                          new_coordinates[2], particle_offset[0],
                          particle_offset[1], particle_offset[2]);
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  ASSERT_EQ(LennardJones::directionalForceWithOffset(
                *p, *q, accounted_particle_distance),
            f.directionalForce(*test, *q));
}

/**
 * tests that the correct force ist calculated for a periodic boundary in the x
 * and y dimension
 */
TEST(PeriodicBoundaryForce, offsetXY1) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({8, 8, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({-1, -1, 1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting(
          ee,
          xhigh);  // xhigh and not yhigh because corners are invalid for yhigh
                   // if both boundaries are periodic to avoid double matching
  constexpr dvec3 expected_offset = {9, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  DVEC3_NEAR(LennardJones::directionalForceWithOffset(
                 *p, *q, accounted_particle_distance),
             f.directionalForce(*test, *q), "not equal to theoretical force",
             1e-8);
}

/**
 * tests that the correct force ist calculated for a periodic boundary in the x
 * and y dimension
 */
TEST(PeriodicBoundaryForce, offsetXY2) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 1},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({8, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 8, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({-1, 10, 1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting(
          ae,
          xhigh);  // xhigh and not yhigh because corners are invalid for yhigh
  // if both boundaries are periodic to avoid double matching
  constexpr dvec3 expected_offset = {9, -9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  DVEC3_NEAR(LennardJones::directionalForceWithOffset(
                 *p, *q, accounted_particle_distance),
             f.directionalForce(*test, *q), "not equal to theoretical force",
             1e-8);
}

/**
 * note that the offset tests do not add any real value over the warp tests,
 * therefore, for the extension to 3D, these will include fewer examples for Z.
 * It only has to be tested that the directionalForceWithOffset accepts a Z
 * offset
 */

/**
 * tests that the correct force ist calculated for a periodic boundary in the x,
 * y and z dimension
 */
TEST(PeriodicBoundaryForce, offsetXYZ) {
  const LinkedCellsContainer container(
      {.domain = {9, 9, 9},
       .cutoff_radius = 3,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
       }});

  const LennardJones f{};

  // distance of 2, test particle is in theoretical position
  const auto p = new Particle({8, 8, 8}, {0, 0, 0}, 1, 5.0, 1);  // in bd
  const auto q = new Particle({1, 1, 1}, {0, 0, 0}, 1, 5.0, 1);  // in bb
  const auto test = new Particle({-1, -1, -1}, {0, 0, 0}, 1, 5.0, 1);

  auto [is_adjacent, new_coordinates, particle_offset] =
      container.reflectiveWarpAroundTesting(
          {3, 3, 3},
          xhigh);  // xhigh and not yhigh or zhigh because corners are invalid
                   // for yhigh and zhigh
  // if all three boundaries are periodic to avoid double matching
  constexpr dvec3 expected_offset = {9, 9, 9};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(particle_offset, expected_offset, "Offset wrong", 1e-5);

  // compare that the force with offset is equal to what is the force in theory
  const dvec3 accounted_particle_distance =
      q->getX() - p->getX() + particle_offset;
  DVEC3_NEAR(LennardJones::directionalForceWithOffset(
                 *p, *q, accounted_particle_distance),
             f.directionalForce(*test, *q), "not equal to theoretical force",
             1e-8);
}

//[[======================== moving test =====================================]]

/**
 * tests that particles can move through the x periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXLeft) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {-2, 0, 0}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {-2.0, 0.0, 0.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {8, 1, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the x periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXRight) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end

  Particle one({8, 1, 1}, {2, 0, 0}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {2.0, 0.0, 0.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the x periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXDiagonal1) {
  LinkedCellsContainer container({.domain = {9, 9, 9},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 8, 1}, {2, 2, 0}, 1, 5.0, 1);
  Particle two({8, 1, 1}, {2, -2, 0}, 1, 5.0, 1);

  std::cout << "Adding Particles" << std::endl;
  container.addParticles({one, two});
  EXPECT_EQ(container.size(), 2) << "Number of Particles is not 2";

  SpdWrapper::get()->info("Starting test");
  std::cout << "Starting Test" << std::endl;

  double delta_t = 1;
  container.singleIterator([this, delta_t](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
    std::cout << "setting p " << p.getId() << " new_x [" << new_x[0] << ","
              << new_x[1] << "," << new_x[2] << "]" << std::endl;
  });

  SpdWrapper::get()->info("Positions Updated");
  std::cout << "Positions Updated" << std::endl;

  container.imposeInvariant();
  SpdWrapper::get()->info("Invariant Imposed");
  std::cout << "Invariant Imposed" << std::endl;

  for (const auto& c : container.getCells()) {
    for (const auto p : c) {
      SpdWrapper::get()->info("p {} is at [{}, {}, {}]", p->getId(),
                              p->getX()[0], p->getX()[1], p->getX()[2]);
      FAIL() << "Particle should have been deleted";
    }
  }
}

/**
 * tests that particles can move through the x periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXDiagonal2) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {-2, -2, 0}, 1, 5.0, 1);
  Particle two({1, 8, 1}, {-2, 2, 0}, 1, 5.0, 1);

  container.addParticles({one, two});
  EXPECT_EQ(container.size(), 2) << "Number of Particles is not 2";

  double delta_t = 1;
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

  for (const auto& c : container.getCells()) {
    for (const auto p : c) {
      SpdWrapper::get()->info("p {} is at [{}, {}, {}]", p->getId(),
                              p->getX()[0], p->getX()[1], p->getX()[2]);
      FAIL() << "Particle should have been deleted";
    }
  }
}

/**
 * tests that particles can move through the y periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveYDiagonal1) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 8, 1}, {2, 2, 0}, 1, 5.0, 1);
  Particle two({8, 1, 1}, {2, -2, 0}, 1, 5.0, 1);

  container.addParticles({one, two});
  EXPECT_EQ(container.size(), 2) << "Number of Particles is not 2";

  double delta_t = 1;
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

  for (const auto& c : container.getCells()) {
    for (const auto p : c) {
      SpdWrapper::get()->info("p {} is at [{}, {}, {}]", p->getId(),
                              p->getX()[0], p->getX()[1], p->getX()[2]);
      FAIL() << "Particle should have been deleted";
    }
  }
}

/**
 * tests that particles can move through the y periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveYDiagonal2) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {-2, -2, 0}, 1, 5.0, 1);
  Particle two({1, 8, 1}, {-2, 2, 0}, 1, 5.0, 1);

  container.addParticles({one, two});
  EXPECT_EQ(container.size(), 2) << "Number of Particles is not 2";

  double delta_t = 1;
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

  for (const auto& c : container.getCells()) {
    for (const auto p : c) {
      SpdWrapper::get()->info("p {} is at [{}, {}, {}]", p->getId(),
                              p->getX()[0], p->getX()[1], p->getX()[2]);
      FAIL() << "Particle should have been deleted";
    }
  }
}

/**
 * tests that particles can move through the y periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveYDown) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {0, -2, 0}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {0, -2.0, 0.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 8, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 2, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the y periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveYUp) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 8, 1}, {0, 2, 0}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {0.0, 2.0, 0.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 2, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the xy periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXYDiagonal1) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {-2, -2, 0}, 1, 5.0, 1);

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getX(), {8, 8, 1}, "position wrong", 1e-5);
    DVEC3_NEAR(p.getV(), {-2, -2, 0}, "velocity wrong", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 2, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the xy periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXYDiagonal2) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 1, 1}, {2, -2, 0}, 1, 5.0, 1);

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getX(), {1, 8, 1}, "position wrong", 1e-5);
    DVEC3_NEAR(p.getV(), {2, -2, 0}, "velocity wrong", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 2, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the xy periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXYDiagonal3) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 8, 1}, {-2, 2, 0}, 1, 5.0, 1);

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getX(), {8, 1, 1}, "position wrong", 1e-5);
    DVEC3_NEAR(p.getV(), {-2, 2, 0}, "velocity wrong", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 2, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the xy periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXYDiagonal4) {
  LinkedCellsContainer container({.domain = {9, 9, 3},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 8, 1}, {2, 2, 0}, 1, 5.0, 1);

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "position wrong", 1e-5);
    DVEC3_NEAR(p.getV(), {2, 2, 0}, "velocity wrong", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 2, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the y periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveZForward) {
  LinkedCellsContainer container({.domain = {9, 9, 9},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 1}, {0, 0, -2}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 2";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {0.0, 0.0, -2.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 1, 8}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 2})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the z periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveZBackward) {
  LinkedCellsContainer container({.domain = {9, 9, 9},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Outflow,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({1, 1, 8}, {0, 0, 2}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {0.0, 0.0, 2.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 2})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}

/**
 * tests that particles can move through the x, y and z periodic boundary
 */
TEST(PeriodicBoundaryMoving, moveXYZDiagonal) {
  LinkedCellsContainer container({.domain = {9, 9, 9},
                                  .cutoff_radius = 3,
                                  .boundary_config = {
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                      LinkedCellsConfig::BoundaryType::Periodic,
                                  }});

  const LennardJones f{};

  // move to other end
  Particle one({8, 8, 8}, {2, 2, 2}, 1, 5.0, 1);  // in bb

  container.addParticles({one});
  EXPECT_EQ(container.size(), 1) << "Number of Particles is not 1";

  double delta_t = 1;
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

  container.singleIterator([this](Particle& p) {
    DVEC3_NEAR(p.getV(), {2.0, 2.0, 2.0}, "wrong velocity", 1e-5);
    DVEC3_NEAR(p.getX(), {1, 1, 1}, "wrong position", 1e-5);
  });

  // test that particle is registered in its true cell
  const auto old_cell =
      container.getCells()[container.cellCoordToIndexTesting({2, 2, 2})];
  const auto new_cell =
      container.getCells()[container.cellCoordToIndexTesting({0, 0, 0})];

  EXPECT_EQ(old_cell.size(), 0);
  EXPECT_EQ(new_cell.size(), 1);
}