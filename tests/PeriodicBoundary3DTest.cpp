//
// Created by maximilian on 22.01.25.
//
#include <gtest/gtest.h>

#include "defs/Particle.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "forces/LennardJones.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"


/**
 * note by the author: this does test ALL combinations for a 3x3 (and therefore 5x5) cube.
 * I spent too much time doing this
 */

/** side view
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

//[[======================== warping tests ===================================]]

/**
 * tests warping in 3D for only x periodic
 * and only the most interesting cells, because else it would be too much
 */
TEST(PeriodicBoundary3D, warpingX) {
  const LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 1,
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

  // row 3
  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ea";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 0}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "eb";

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ec";

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 2}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ed";

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ee";

  //----------------------------------------//
  // row2
  // da
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "da";

  // db
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 0}, xhigh);
  constexpr dvec3 expected_offset11 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "db";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 0});
  DVEC3_NEAR(offset, expected_offset11, "Offset wrong", 1e-5);

  // dc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 1}, xhigh);
  constexpr dvec3 expected_offset12 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "dc";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 1});
  DVEC3_NEAR(offset, expected_offset12, "Offset wrong", 1e-5);

  // dd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 2}, xhigh);
  constexpr dvec3 expected_offset13 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "dd";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 2});
  DVEC3_NEAR(offset, expected_offset13, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "de";

  //----------------------------------------//
  // row1
  // ca
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ca";

  // cb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 0}, xhigh);
  constexpr dvec3 expected_offset21 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "cb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 0});
  DVEC3_NEAR(offset, expected_offset21, "Offset wrong", 1e-5);

  // cc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 1}, xhigh);
  constexpr dvec3 expected_offset22 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "cc";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 1});
  DVEC3_NEAR(offset, expected_offset22, "Offset wrong", 1e-5);

  // cd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 2}, xhigh);
  constexpr dvec3 expected_offset23 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "cd";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 2});
  DVEC3_NEAR(offset, expected_offset23, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ce";

  //----------------------------------------//
  // row0
  // ba
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ba";

  // bb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 0}, xhigh);
  constexpr dvec3 expected_offset31 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "bb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset31, "Offset wrong", 1e-5);

  // bc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 1}, xhigh);
  constexpr dvec3 expected_offset32 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "bc";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 1});
  DVEC3_NEAR(offset, expected_offset32, "Offset wrong", 1e-5);

  // bd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 2}, xhigh);
  constexpr dvec3 expected_offset33 = {3, 0, 0};
  EXPECT_EQ(is_adjacent, true) << "bd";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 2});
  DVEC3_NEAR(offset, expected_offset33, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "be";

  //----------------------------------------//
  // row-1
  // aa
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "aa";

  // ab
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 0}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ab";

  // ac
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ac";

  // ad
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 2}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ad";

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ae";
}

/**
 * tests warping in 3D for only y periodic
 * and only the most interesting cells, because else it would be too much
 */
TEST(PeriodicBoundary3D, warpingY) {
  const LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 1,
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

  // row 3
  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 3, 3}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ea";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 3}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "eb";

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 3}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ec";

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 3}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ed";

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 3}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ee";

  //----------------------------------------//
  // row2
  // da
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 3, 2}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "da";

  // db
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 2}, yhigh);
  constexpr dvec3 expected_offset11 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "db";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 2});
  DVEC3_NEAR(offset, expected_offset11, "Offset wrong", 1e-5);

  // dc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 2}, yhigh);
  constexpr dvec3 expected_offset12 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "dc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 2});
  DVEC3_NEAR(offset, expected_offset12, "Offset wrong", 1e-5);

  // dd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 2}, yhigh);
  constexpr dvec3 expected_offset13 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "dd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 2});
  DVEC3_NEAR(offset, expected_offset13, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 2}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "de";

  //----------------------------------------//
  // row1
  // ca
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 3, 1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ca";

  // cb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 1}, yhigh);
  constexpr dvec3 expected_offset21 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 1});
  DVEC3_NEAR(offset, expected_offset21, "Offset wrong", 1e-5);

  // cc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 1}, yhigh);
  constexpr dvec3 expected_offset22 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 1});
  DVEC3_NEAR(offset, expected_offset22, "Offset wrong", 1e-5);

  // cd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 1}, yhigh);
  constexpr dvec3 expected_offset23 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 1});
  DVEC3_NEAR(offset, expected_offset23, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ce";

  //----------------------------------------//
  // row0
  // ba
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 3, 0}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ba";

  // bb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 0}, yhigh);
  constexpr dvec3 expected_offset31 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset31, "Offset wrong", 1e-5);

  // bc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 0}, yhigh);
  constexpr dvec3 expected_offset32 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 0});
  DVEC3_NEAR(offset, expected_offset32, "Offset wrong", 1e-5);

  // bd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 0}, yhigh);
  constexpr dvec3 expected_offset33 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 0});
  DVEC3_NEAR(offset, expected_offset33, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 0}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "be";

  //----------------------------------------//
  // row-1
  // aa
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 3, -1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "aa";

  // ab
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, -1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ab";

  // ac
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, -1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ac";

  // ad
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, -1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ad";

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, -1}, yhigh);
  EXPECT_EQ(is_adjacent, false) << "ae";
}