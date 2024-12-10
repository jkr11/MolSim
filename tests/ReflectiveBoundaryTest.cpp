//
// Created by maximilian on 10.12.24.
//
#include <gtest/gtest.h>

#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
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

/**
 * tests that warping is correct for x cells in 2D for a 3x3 cell grid
 * with both axis Periodic
 */
TEST(ReflectiveBoundary, warpingBothPeriodicForX) {
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
      container.reflective_warp_around_testing(ee, xhigh);
  constexpr dvec3 expected_offset1 = {9, 9, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset1, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(de, xhigh);
  constexpr dvec3 expected_offset2 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ce, xhigh);
  constexpr dvec3 expected_offset3 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, cb);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(be, xhigh);
  constexpr dvec3 expected_offset4 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ae, xhigh);
  constexpr dvec3 expected_offset5 = {9, -9, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset5, "Offset wrong", 1e-5);
}

/**
 * tests that warping is correct for y cells in 2D for a 3x3 cell grid
 * with both axis Periodic
 */
TEST(ReflectiveBoundary, warpingBothPeriodicForY) {
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
      container.reflective_warp_around_testing(ea, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency1";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(eb, yhigh);
  constexpr dvec3 expected_offset2 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency2";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ec, yhigh);
  constexpr dvec3 expected_offset3 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency3";
  EXPECT_IVEC3_EQ(new_coordinates, bc);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ed, yhigh);
  constexpr dvec3 expected_offset4 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency4";
  EXPECT_IVEC3_EQ(new_coordinates, bd);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ee, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency5";
}

/**
 * tests that warping is correct for x cells in 2D for a 3x3 cell grid
 * with x axis Periodic
 */
TEST(ReflectiveBoundary, warpingXPeriodicForX) {
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
      container.reflective_warp_around_testing(ee, xhigh);
  EXPECT_EQ(is_adjacent, false);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(de, xhigh);
  constexpr dvec3 expected_offset2 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, db);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ce, xhigh);
  constexpr dvec3 expected_offset3 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, cb);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(be, xhigh);
  constexpr dvec3 expected_offset4 = {9, 0, 0};
  EXPECT_EQ(is_adjacent, true);
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ae, xhigh);
  EXPECT_EQ(is_adjacent, false);
}

/**
 * tests that warping is correct for y cells in 2D for a 3x3 cell grid
 * with y axis Periodic
 */
TEST(ReflectiveBoundary, warpingYPeriodicForY) {
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
      container.reflective_warp_around_testing(ea, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency1";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(eb, yhigh);
  constexpr dvec3 expected_offset2 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency2";
  EXPECT_IVEC3_EQ(new_coordinates, bb);
  DVEC3_NEAR(offset, expected_offset2, "Offset wrong", 1e-5);

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ec, yhigh);
  constexpr dvec3 expected_offset3 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency3";
  EXPECT_IVEC3_EQ(new_coordinates, bc);
  DVEC3_NEAR(offset, expected_offset3, "Offset wrong", 1e-5);

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ed, yhigh);
  constexpr dvec3 expected_offset4 = {0, 9, 0};
  EXPECT_EQ(is_adjacent, true) << "wrong adjacency4";
  EXPECT_IVEC3_EQ(new_coordinates, bd);
  DVEC3_NEAR(offset, expected_offset4, "Offset wrong", 1e-5);

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing(ee, yhigh);
  EXPECT_EQ(is_adjacent, false) << "wrong adjacency5";
}


