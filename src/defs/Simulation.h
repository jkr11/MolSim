//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <variant>

#include "defs/types.h"
#include "utils/SpdWrapper.h"

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

struct SingularGravityConfig {
  double g{};
};

struct ThermostatConfig {
  double T_init{};
  double T_target{};
  double deltaT{};
  int n_thermostat{};
};
// TODO: apparently we cant nest these for access in XMLReader
/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  double t_end;
  double delta_t;
  enum ForceType { LennardJones, Gravity } force_type;
  enum SingularForceType { SingularGravity } singular_force_type;
  ThermostatConfig thermostat_config;
  // TODO: remove this vvvvv
  std::variant<SingularGravityConfig> singular_force_data;
  std::variant<LinkedCellsConfig, DirectSumConfig> container_data;
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

  logger->info("Force Type: {}", args.force_type == Arguments::LennardJones
                                     ? "Lennard-Jones"
                                     : "Gravity");

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
