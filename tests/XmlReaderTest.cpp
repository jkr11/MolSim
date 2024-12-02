#include <gtest/gtest.h>

#include <fstream>

#include "../src/io/file/in/xml/XmlReader.h"
#include "../src/io/file/in/xml/input.cxx"  // It wants this idk why
#include "../src/io/file/in/xml/input.hxx"
#include "spdlog/fmt/bundled/os.h"
#include "testUtil.h"

Arguments arguments = {
    .t_end = 5,
    .delta_t = 0.0002,
    .force_type = Arguments::LennardJones,
    .container_data =
        LinkedCellsConfig{.domain = {100, 100, 100},
                          .cutoff_radius = 3.0,
                          .boundary_config =
                              {
                                  .x_high = LinkedCellsConfig::Outflow,
                                  .x_low = LinkedCellsConfig::Outflow,
                                  .y_high = LinkedCellsConfig::Outflow,
                                  .y_low = LinkedCellsConfig::Outflow,
                                  .z_high = LinkedCellsConfig::Outflow,
                                  .z_low = LinkedCellsConfig::Outflow,
                              }},
};

TEST(XmlReader, failOnNonXml) {
  if (std::ofstream file("testFile.txt"); file.is_open()) {
    file << "Hello, World!" << std::endl;
    file.close();
  }
  std::vector<Particle> particles;
  Arguments args = arguments;
  EXPECT_ANY_THROW(XmlReader::read(particles, "testFile.txt", args));
}

TEST(XmlReader, testCuboid) {
  std::vector<Particle> particles;
  XmlReader::read(particles, "tests/test_cuboid.xml", arguments);

  EXPECT_EQ(particles.size(), 20);
  EXPECT_EQ(arguments.t_end, 15);
  EXPECT_EQ(arguments.delta_t, 0.015);
  EXPECT_EQ(arguments.force_type, Arguments::LennardJones);
  EXPECT_TRUE(
      std::holds_alternative<DirectSumConfig>(arguments.container_data));
}

TEST(XmlReader, testCuboidSpheroidLinkedCells) {
  std::vector<Particle> particles;
  XmlReader::read(particles, "tests/test_cuboid_spheroid.xml", arguments);
  EXPECT_EQ(particles.size(), 1515);
  EXPECT_EQ(arguments.t_end, 15);
  EXPECT_EQ(arguments.delta_t, 0.015);
  EXPECT_EQ(arguments.force_type, Arguments::LennardJones);
  EXPECT_TRUE(
      std::holds_alternative<LinkedCellsConfig>(arguments.container_data));
  const auto &config = std::get<LinkedCellsConfig>(arguments.container_data);
  constexpr auto comp = ivec3{400, 400, 1};
  EXPECT_IVEC3_EQ(config.domain, comp);
  EXPECT_EQ(config.cutoff_radius, 3.0);
  EXPECT_EQ(config.boundary_config.x_high, LinkedCellsConfig::Outflow);
  EXPECT_EQ(config.boundary_config.x_low, LinkedCellsConfig::Outflow);
  EXPECT_EQ(config.boundary_config.y_high, LinkedCellsConfig::Outflow);
  EXPECT_EQ(config.boundary_config.y_low, LinkedCellsConfig::Outflow);
  EXPECT_EQ(config.boundary_config.z_high, LinkedCellsConfig::Outflow);
  EXPECT_EQ(config.boundary_config.z_low, LinkedCellsConfig::Outflow);
}