#include <gtest/gtest.h>

#include <fstream>

#include "../src/io/file/in/xml/XmlReader.h"
#include "io/file/in/xml/input.cxx"
#include "io/file/in/xml/input.hxx"
#include "io/file/out/checkpoint-schema.cxx"
#include "io/file/out/checkpoint-schema.hxx"
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

/**
 *  @brief tests if xmlReader fails for non xml files
 */
TEST(XmlReader, failOnNonXml) {
  if (std::ofstream file("testFile.txt"); file.is_open()) {
    file << "Hello, World!" << std::endl;
    file.close();
  }
  std::vector<Particle> particles;
  Arguments args = arguments;
  EXPECT_ANY_THROW(XmlReader::read(particles, "testFile.txt", args));
}

/**
 *  @brief checks if cuboids and DirectSum containers are read correctly
 *  this test will only pass if it is executed via our build script due to
 * hardcoded file paths
 */
TEST(XmlReader, testCuboid) {
  std::vector<Particle> particles;
  XmlReader::read(particles, "../../tests/test_cuboid.xml", arguments);

  EXPECT_EQ(particles.size(), 20);
  EXPECT_EQ(arguments.t_end, 15);
  EXPECT_EQ(arguments.delta_t, 0.015);
  EXPECT_EQ(arguments.force_type, Arguments::LennardJones);
  EXPECT_TRUE(
      std::holds_alternative<DirectSumConfig>(arguments.container_data));
}

/**
 * @brief checks if both cuboids and spheroids are correct, LinkedCellsConfig
 * this test will only pass if it is executed via our build script due to
 * hardcoded file paths
 */
TEST(XmlReader, testCuboidSpheroidLinkedCells) {
  std::vector<Particle> particles;
  XmlReader::read(particles, "../../tests/test_cuboid_spheroid.xml", arguments);
  EXPECT_EQ(particles.size(), 1257);
  EXPECT_EQ(arguments.t_end, 15);
  EXPECT_EQ(arguments.delta_t, 0.015);
  EXPECT_EQ(arguments.force_type, Arguments::LennardJones);
  EXPECT_TRUE(
      std::holds_alternative<LinkedCellsConfig>(arguments.container_data));
  const auto& config = std::get<LinkedCellsConfig>(arguments.container_data);
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

namespace fs = std::filesystem;

[[nodiscard]] bool isXMLFile(const fs::path& filePath) noexcept {
  return filePath.extension() == ".xml";
}

void processXMLFilesInInput(std::vector<fs::path>& paths) {
  const std::string inputDir = "../../input";

  try {
    if (!fs::exists(inputDir)) {
      std::cerr << "Input directory does not exist: " << inputDir << std::endl;
      return;
    }
    for (const auto& entry : fs::directory_iterator(inputDir)) {
      if (entry.is_regular_file() && isXMLFile(entry.path()) &&
          !fs::is_directory(entry.path())) {
        paths.emplace_back(entry.path());
      }
    }
  } catch (const fs::filesystem_error& e) {
    std::cerr << "Filesystem error: " << e.what() << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}
TEST(XmlReader, all_inputs_no_error) {
  std::vector<fs::path> paths;
  processXMLFilesInInput(paths);

  for (const auto& path : paths) {
    if (path != "../../input/week43.xml" &&
        path != "../../input/week43periodic.xml" &&
        path != "../../input/checkpoint_test.xml") {  // because this requires
                                                      // checkpoint
      SpdWrapper::get()->info("Path {}", path.c_str());
      // that isnt written yet
      std::vector<Particle> particles;
      EXPECT_NO_FATAL_FAILURE(XmlReader::read(particles, path, arguments));
    }
  }
}