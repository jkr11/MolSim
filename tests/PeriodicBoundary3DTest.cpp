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
 * note by the author: this does test ALL combinations for a 3x3 (and therefore
 * 5x5) cube. I spent too much time doing this
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
/// -------------------- [[ 1 boundary periodic ]] -------------------- ///

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

/**
 * tests warping in 3D for only z periodic
 * and only the most interesting cells, because else it would be too much
 */
TEST(PeriodicBoundary3D, warpingZ) {
  const LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // row 3
  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ea";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "eb";

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ec";

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ed";

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ee";

  //----------------------------------------//
  // row2
  // da
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "da";

  // db
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 2, 3}, zhigh);
  constexpr dvec3 expected_offset11 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "db";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 0});
  DVEC3_NEAR(offset, expected_offset11, "Offset wrong", 1e-5);

  // dc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 2, 3}, zhigh);
  constexpr dvec3 expected_offset12 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "dc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 2, 0});
  DVEC3_NEAR(offset, expected_offset12, "Offset wrong", 1e-5);

  // dd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 2, 3}, zhigh);
  constexpr dvec3 expected_offset13 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "dd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 2, 0});
  DVEC3_NEAR(offset, expected_offset13, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "de";

  //----------------------------------------//
  // row1
  // ca
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ca";

  // cb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 1, 3}, zhigh);
  constexpr dvec3 expected_offset21 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 0});
  DVEC3_NEAR(offset, expected_offset21, "Offset wrong", 1e-5);

  // cc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 1, 3}, zhigh);
  constexpr dvec3 expected_offset22 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 1, 0});
  DVEC3_NEAR(offset, expected_offset22, "Offset wrong", 1e-5);

  // cd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 1, 3}, zhigh);
  constexpr dvec3 expected_offset23 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 1, 0});
  DVEC3_NEAR(offset, expected_offset23, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ce";

  //----------------------------------------//
  // row0
  // ba
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 0, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ba";

  // bb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 0, 3}, zhigh);
  constexpr dvec3 expected_offset31 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset31, "Offset wrong", 1e-5);

  // bc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 0, 3}, zhigh);
  constexpr dvec3 expected_offset32 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 0});
  DVEC3_NEAR(offset, expected_offset32, "Offset wrong", 1e-5);

  // bd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 0, 3}, zhigh);
  constexpr dvec3 expected_offset33 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 0});
  DVEC3_NEAR(offset, expected_offset33, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "be";

  //----------------------------------------//
  // row-1
  // aa
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "aa";

  // ab
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ab";

  // ac
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ac";

  // ad
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ad";

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ae";
}

/// -------------------- [[ 2 boundaries periodic ]] -------------------- ///
/**
 * tests warping in 3D for only x and y periodic
 * and only the most interesting cells, because else it would be too much
 */
TEST(PeriodicBoundary3D, warpingXY) {
  const LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 1,
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

  // on x
  // row 3
  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, -1}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ea";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 0}, xhigh);
  constexpr dvec3 expected_offset01 = {3, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "eb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset01, "Offset wrong", 1e-5);

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 1}, xhigh);
  constexpr dvec3 expected_offset02 = {3, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "ec";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 1});
  DVEC3_NEAR(offset, expected_offset02, "Offset wrong", 1e-5);

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 2}, xhigh);
  constexpr dvec3 expected_offset03 = {3, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "ed";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 2});
  DVEC3_NEAR(offset, expected_offset03, "Offset wrong", 1e-5);

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
  constexpr dvec3 expected_offset41 = {3, -3, 0};
  EXPECT_EQ(is_adjacent, true) << "ab";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 0});
  DVEC3_NEAR(offset, expected_offset41, "Offset wrong", 1e-5);

  // ac
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 1}, xhigh);
  constexpr dvec3 expected_offset42 = {3, -3, 0};
  EXPECT_EQ(is_adjacent, true) << "ac";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 1});
  DVEC3_NEAR(offset, expected_offset42, "Offset wrong", 1e-5);

  // ad
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 2}, xhigh);
  constexpr dvec3 expected_offset43 = {3, -3, 0};
  EXPECT_EQ(is_adjacent, true) << "ad";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 2});
  DVEC3_NEAR(offset, expected_offset43, "Offset wrong", 1e-5);

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 3}, xhigh);
  EXPECT_EQ(is_adjacent, false) << "ae";

  // ---------------------------------------------------------------------------
  // on y
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
  constexpr dvec3 expected_offset111 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "db";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 2});
  DVEC3_NEAR(offset, expected_offset111, "Offset wrong", 1e-5);

  // dc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 2}, yhigh);
  constexpr dvec3 expected_offset112 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "dc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 2});
  DVEC3_NEAR(offset, expected_offset112, "Offset wrong", 1e-5);

  // dd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 2}, yhigh);
  constexpr dvec3 expected_offset113 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "dd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 2});
  DVEC3_NEAR(offset, expected_offset113, "Offset wrong", 1e-5);

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
  constexpr dvec3 expected_offset121 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 1});
  DVEC3_NEAR(offset, expected_offset121, "Offset wrong", 1e-5);

  // cc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 1}, yhigh);
  constexpr dvec3 expected_offset122 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 1});
  DVEC3_NEAR(offset, expected_offset122, "Offset wrong", 1e-5);

  // cd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 1}, yhigh);
  constexpr dvec3 expected_offset123 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "cd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 1});
  DVEC3_NEAR(offset, expected_offset123, "Offset wrong", 1e-5);

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
  constexpr dvec3 expected_offset131 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset131, "Offset wrong", 1e-5);

  // bc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 0}, yhigh);
  constexpr dvec3 expected_offset132 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 0});
  DVEC3_NEAR(offset, expected_offset132, "Offset wrong", 1e-5);

  // bd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 0}, yhigh);
  constexpr dvec3 expected_offset133 = {0, 3, 0};
  EXPECT_EQ(is_adjacent, true) << "bd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 0});
  DVEC3_NEAR(offset, expected_offset133, "Offset wrong", 1e-5);

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

/**
 * tests warping in 3D for only x and z periodic
 * and only the most interesting cells, because else it would be too much
 */
TEST(PeriodicBoundary3D, warpingXZ) {
  const LinkedCellsContainer container(
      {.domain = {3, 3, 3},
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Periodic,
           LinkedCellsConfig::BoundaryType::Periodic,
       }});
  bool is_adjacent;
  ivec3 new_coordinates;
  dvec3 offset;

  // on x
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
  constexpr dvec3 expected_offset10 = {3, 0, -3};
  EXPECT_EQ(is_adjacent, true) << "da";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 2});
  DVEC3_NEAR(offset, expected_offset10, "Offset wrong", 1e-5);

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
  constexpr dvec3 expected_offset14 = {3, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "de";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 0});
  DVEC3_NEAR(offset, expected_offset14, "Offset wrong", 1e-5);

  //----------------------------------------//
  // row1
  // ca
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, -1}, xhigh);
  constexpr dvec3 expected_offset20 = {3, 0, -3};
  EXPECT_EQ(is_adjacent, true) << "ca";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 2});
  DVEC3_NEAR(offset, expected_offset20, "Offset wrong", 1e-5);

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
  constexpr dvec3 expected_offset24 = {3, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "ce";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 0});
  DVEC3_NEAR(offset, expected_offset24, "Offset wrong", 1e-5);

  //----------------------------------------//
  // row0
  // ba
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, -1}, xhigh);
  constexpr dvec3 expected_offset30 = {3, 0, -3};
  EXPECT_EQ(is_adjacent, true) << "ba";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 2});
  DVEC3_NEAR(offset, expected_offset30, "Offset wrong", 1e-5);
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
  constexpr dvec3 expected_offset34 = {3, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "be";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset34, "Offset wrong", 1e-5);

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

  // ---------------------------------------------------------------------------
  // on z
  // row 3
  // ea
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ea";

  // eb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "eb";

  // ec
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ec";

  // ed
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ed";

  // ee
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 3, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ee";

  //----------------------------------------//
  // row2
  // da
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "da";

  // db
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 2, 3}, zhigh);
  constexpr dvec3 expected_offset111 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "db";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 2, 0});
  DVEC3_NEAR(offset, expected_offset111, "Offset wrong", 1e-5);

  // dc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 2, 3}, zhigh);
  constexpr dvec3 expected_offset112 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "dc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 2, 0});
  DVEC3_NEAR(offset, expected_offset112, "Offset wrong", 1e-5);

  // dd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 2, 3}, zhigh);
  constexpr dvec3 expected_offset113 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "dd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 2, 0});
  DVEC3_NEAR(offset, expected_offset113, "Offset wrong", 1e-5);

  // de
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 2, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "de";

  //----------------------------------------//
  // row1
  // ca
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ca";

  // cb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 1, 3}, zhigh);
  constexpr dvec3 expected_offset121 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 1, 0});
  DVEC3_NEAR(offset, expected_offset121, "Offset wrong", 1e-5);

  // cc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 1, 3}, zhigh);
  constexpr dvec3 expected_offset122 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 1, 0});
  DVEC3_NEAR(offset, expected_offset122, "Offset wrong", 1e-5);

  // cd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 1, 3}, zhigh);
  constexpr dvec3 expected_offset123 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "cd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 1, 0});
  DVEC3_NEAR(offset, expected_offset123, "Offset wrong", 1e-5);

  // ce
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ce";

  //----------------------------------------//
  // row0
  // ba
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, 0, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ba";

  // bb
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, 0, 3}, zhigh);
  constexpr dvec3 expected_offset131 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bb";
  EXPECT_IVEC3_EQ(new_coordinates, {0, 0, 0});
  DVEC3_NEAR(offset, expected_offset131, "Offset wrong", 1e-5);

  // bc
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, 0, 3}, zhigh);
  constexpr dvec3 expected_offset132 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bc";
  EXPECT_IVEC3_EQ(new_coordinates, {1, 0, 0});
  DVEC3_NEAR(offset, expected_offset132, "Offset wrong", 1e-5);

  // bd
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, 0, 3}, zhigh);
  constexpr dvec3 expected_offset133 = {0, 0, 3};
  EXPECT_EQ(is_adjacent, true) << "bd";
  EXPECT_IVEC3_EQ(new_coordinates, {2, 0, 0});
  DVEC3_NEAR(offset, expected_offset133, "Offset wrong", 1e-5);

  // be
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, 0, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "be";

  //----------------------------------------//
  // row-1
  // aa
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({-1, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "aa";

  // ab
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({0, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ab";

  // ac
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({1, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ac";

  // ad
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({2, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ad";

  // ae
  std::tie(is_adjacent, new_coordinates, offset) =
      container.reflective_warp_around_testing({3, -1, 3}, zhigh);
  EXPECT_EQ(is_adjacent, false) << "ae";
}