//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <variant>

#include "defs/types.h"
#include "forces/SingularForce.h"
#include "utils/SpdWrapper.h"

struct IndexForceConfig {
  std::vector<ivec3> indeces{};
  std::vector<int> ids{};
  double end_time{};
  dvec3 force_values{};
};

/**
 * @brief holds the specification for the LinkedCellsContainer
 */
struct LinkedCellsConfig {
  ivec3 domain;  // size of the dimensions
  double cutoff_radius;
  enum BoundaryType { Outflow, Reflective, Periodic } boundary_type;
  struct BoundaryConfig {
    BoundaryType x_high;
    BoundaryType x_low;
    BoundaryType y_high;
    BoundaryType y_low;
    BoundaryType z_high;
    BoundaryType z_low;
  } boundary_config;
};

/**
 * @brief holds the specification for the DirectSumContainer
 * @note this is empty and probably will always be, but its nice for
 * std::variant
 */
struct DirectSumConfig {};

/**
 * @brief holds the specification for the SingularGravity force
 */
struct SingularGravityConfig {
  double g{};
  int axis{};
};



/**
 * @brief holds the harmonic force parameters
 */
struct HarmonicForceConfig {
  double r_0{};
  double k{};
};

/**
 * @brief holds the LennardJones force parameters
 */
struct LennardJonesConfig {};

struct TruncatedLennardJonesConfig {};

/**
 * @brief holds the interactive gravity parameters
 */
struct GravityConfig {};

/**
 * @brief holds instance data for Thermostat
 */
struct ThermostatConfig {
  double t_init{};
  double t_target{};
  double delta_t{};
  int n_thermostat{};
  bool use_relative{};
  bool use_thermal_motion{};
};

struct StatisticsConfig {
  bool calc_stats{};
  int x_bins{};
  int y_bins{};
  int output_interval{};
  std::string velocity_output_location{};
  std::string density_output_location{};
};

struct SphereoidGeneratorConfig {
  dvec3 origin{};
  const int radius{};
  double h{};
  double m{};
  const dvec3 initial_velocity{};
  double epsilon{};
  double sigma{};
  const int type{};
  double mv{};
  const bool two_d{};
  SphereoidGeneratorConfig() = default;


};

struct CuboidGeneratorConfig {
  dvec3 corner{};
  ivec3 dimensions{};
  double h{};
  double m{};
  const dvec3 initial_velocity{};
  double mv{};
  double epsilon{};
  double sigma{};
  const int type{};
  const bool two_d{};
  CuboidGeneratorConfig() = default;

};

struct MembraneGeneratorConfig {
  dvec3 corner{};
  ivec3 dimensions{};
  double h{};
  double m{};
  const dvec3 initial_velocity{};
  double mv{};
  double epsilon{};
  double sigma{};
  const int type{};
  const bool two_d{};
  std::vector<int> ids{};
  std::vector<ivec3> indeces{};
};

/**
 * @brief struct to hold command line arguments
 */
using SingularForceTypes =
    std::variant<SingularGravityConfig, HarmonicForceConfig>;
using InteractiveForceTypes = std::variant<LennardJonesConfig, GravityConfig,
                                           TruncatedLennardJonesConfig>;
struct Arguments {
  double t_end;
  double delta_t;
  ThermostatConfig thermostat_config;
  bool use_thermostat;
  std::variant<LinkedCellsConfig, DirectSumConfig> container_data;
  std::vector<SingularForceTypes> singular_force_types;
  std::vector<InteractiveForceTypes> interactive_force_types;
  std::vector<IndexForceConfig> index_force_configs;
  StatisticsConfig statistics_config;
  SphereoidGeneratorConfig spheroid_generator_config;
  CuboidGeneratorConfig cuboid_generator_config;
  MembraneGeneratorConfig membrane_generator_config;
};

/**
 * @brief converts a complex type to a printable statement for any of the six
 * boundaries and their 3 possible types
 * @param boundary_name name of the boundary that will be printed
 * @param type type of the boundary passed from the struct
 */
inline void printBoundary(const std::string& boundary_name,
                          const LinkedCellsConfig::BoundaryType type) {
  const auto logger = SpdWrapper::get();

  std::string type_str;
  switch (type) {
    case LinkedCellsConfig::Outflow:
      type_str = "Outflow";
      break;
    case LinkedCellsConfig::Reflective:
      type_str = "Reflective";
      break;
    case LinkedCellsConfig::Periodic:
      type_str = "Periodic";
      break;
  }

  logger->info("{} Boundary: {}", boundary_name, type_str);
}

/**
 * @brief prints the configuration data of our simulation, specific to
 * especially the used container
 * @param args the struct for which the config is to be printed
 */
inline void printConfiguration(const Arguments& args) {
  const auto logger = SpdWrapper::get();

  logger->info("Simulation Configuration:");
  logger->info("============================");
  logger->info("t_end: {}", args.t_end);
  logger->info("delta_t: {}", args.delta_t);
  logger->info("Singular Forces:");
  if (args.use_thermostat) {
    logger->info("Thermostat: T_init {}", args.thermostat_config.t_init);
    logger->info("--- T_target: {}", args.thermostat_config.t_target);
    logger->info("--- deltaT: {}", args.thermostat_config.delta_t);
  }

  if (std::holds_alternative<LinkedCellsConfig>(args.container_data)) {
    logger->info("Container Type: Linked Cells");

    const auto& [domain, cutoff_radius, boundary_type, boundary_config] =
        std::get<LinkedCellsConfig>(args.container_data);
    logger->info("-- Domain: ({}, {}, {})", domain[0], domain[1], domain[2]);
    logger->info("-- Cutoff Radius: {}", cutoff_radius);
    logger->info("Boundary Configuration:");
    logger->info("------------------------");

    const auto& [x_high, x_low, y_high, y_low, z_high, z_low] = boundary_config;
    printBoundary("x_high", x_high);
    printBoundary("x_low", x_low);
    printBoundary("y_high", y_high);
    printBoundary("y_low", y_low);
    printBoundary("z_high", z_high);
    printBoundary("z_low", z_low);

  } else {
    logger->info("Container Type: Direct Sum");
  }

  logger->info("============================");
}

#endif  // SIMULATION_H
