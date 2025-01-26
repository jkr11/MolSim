//
// Created by maximilian on 22.01.25.
//
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "defs/Particle.h"
#include "defs/Thermostat.h"
#include "defs/containers/LinkedCellsContainer.h"
#include "testUtil.h"
#include "utils/ArrayUtils.h"
#include "utils/Statistics.h"

/**
 * tests that wall particles are excluded from the thermostat
 */
TEST(Statistics, WallThermostat) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  constexpr ThermostatConfig config = {
    .t_init = 0,
    .t_target = 1,
    .delta_t = 1,
    .n_thermostat = 1,
    .use_relative = false,
    .use_thermal_motion = true,
  };

  Thermostat thermostat(config);

  const Particle wall({2, 1, 1}, {1, 1, 1}, 1, 1, 1, -1);
  const Particle p({1, 2, 1}, {2, 2, 2}, 1, 1, 1, 1);
  const Particle q({1, 1, 2}, {3, 3, 3}, 1, 1, 1, 2);
  container.addParticle(wall);
  container.addParticle(p);
  container.addParticle(q);

  EXPECT_EQ(container.getParticleCount(), 3) << "Particle Count wrong";
  EXPECT_EQ(container.getSpecialParticleCount(), 1) << "Special Particle Count wrong";
  EXPECT_EQ(container.size(), 3) << "Number of Particles is not 0";

  DVEC3_NEAR(Thermostat::getAverageVelocity(container), {2.5, 2.5, 2.5}, "average velocity wrong", 1e-5);
}

/**
 * tests that the thermal temperature ignores walls
 */
TEST(Statistics, ThermalTemperature) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  constexpr ThermostatConfig config = {
    .t_init = 0,
    .t_target = 1,
    .delta_t = 1,
    .n_thermostat = 1,
    .use_relative = false,
    .use_thermal_motion = true,
  };

  Thermostat thermostat(config);

  const Particle wall({2, 1, 1}, {1, 1, 1}, 1, 1, 1, -1);
  const Particle p({1, 2, 1}, {2, 2, 2}, 1, 1, 1, 1);
  const Particle q({1, 1, 2}, {3, 3, 3}, 2, 1, 1, 2);
  container.addParticle(wall);
  container.addParticle(p);
  container.addParticle(q);

  EXPECT_EQ(container.getParticleCount(), 3) << "Particle Count wrong";
  EXPECT_EQ(container.getSpecialParticleCount(), 1) << "Special Particle Count wrong";
  EXPECT_EQ(container.size(), 3) << "Number of Particles is not 0";

  dvec3 average_velocity = Thermostat::getAverageVelocity(container);
  DVEC3_NEAR(average_velocity, {2.5, 2.5, 2.5}, "average velocity wrong", 1e-5);

  EXPECT_EQ(thermostat.getThermalTemperature(container, average_velocity), 0.375);
}

/**
 * tests that the statistics are correct
 * tested with ybins = 1, because this is the official bin type
 */
TEST(Statistics, bins) {
  LinkedCellsContainer container(
      {.domain = {3, 3, 3},  // better safe than sorry
       .cutoff_radius = 1,
       .boundary_config = {
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
           LinkedCellsConfig::BoundaryType::Outflow,
       }});

  std::filesystem::path tempDir = std::filesystem::temp_directory_path();
  std::string densityFile = (tempDir / "density.txt").string();
  std::string velocityFile = (tempDir / "velocity.txt").string();

  Statistics statistics(3, 1, container, densityFile, velocityFile);

  const Particle wall({0.5, 0.5, 0.5}, {1, 1, 1}, 1, 1, 1, -1);
  const Particle p({0.5, 0.5, 0.5}, {2, 2, 2}, 1, 1, 1, 1);
  const Particle q({0.5, 0.5, 0.5}, {3, 3, 1}, 2, 1, 1, 2);
  const Particle r({2.5, 2.5, 0.5}, {4, 3, 3}, 2, 1, 1, 2);

  container.addParticle(wall);
  container.addParticle(p);
  container.addParticle(q);
  container.addParticle(r);

  EXPECT_EQ(container.getParticleCount(), 4) << "Particle Count wrong";

  statistics.writeStatistics(1);
  statistics.closeFiles();

  std::ifstream dfile(densityFile);
  std::ifstream vfile(velocityFile);

  if (!dfile || !vfile) {
    FAIL() << "Could not open at least one of the files";
  }

  std::string line;
  std::getline(dfile, line);
  ASSERT_EQ(line, "1,0.222222,0.000000,0.111111");

 std::getline(vfile, line);
  ASSERT_EQ(line, "1,2.500000 2.500000 1.500000,0.0 0.0 0.0,4.000000 3.000000 3.000000");
}