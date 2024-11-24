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
  ivec3 domain;
  double cutoff_radius;
  enum BoundaryType { Outflow, Reflective, Periodic } boundary_type;
  struct BoundaryConfig {
    BoundaryType north;
    BoundaryType south;
    BoundaryType east;
    BoundaryType west;
    BoundaryType up;
    BoundaryType down;
  } boundary_config;
};

/**
 * @brief holds the specification for the DirectSumContainer
 */
struct DirectSumConfig {
  // this doesnt need any other info.
};

struct SingularGravityConfig {
  double g{};
};
// TODO: apparently we cant nest these for access in XMLReader
/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  double t_end;
  double delta_t;
  // double output_time_step_size; // TODO: should this go in files?
  // TODO: i think log level can be removed too
  std::string log_level;
  enum ForceType { LennardJones, Gravity } force_type;
  enum SingularForceType { SingularGravity } singular_force_type;
  // TODO: remove this vvvvv
  enum ContainerType { LinkedCells, DirectSum } container_type;
  std::variant<SingularGravityConfig> singular_force_data;
  std::variant<LinkedCellsConfig, DirectSumConfig> container_data;
  void printConfiguration() const {
    const auto logger = SpdWrapper::get();

    logger->info("Simulation Configuration:");
    logger->info("============================");
    logger->info("t_end: {}", t_end);
    logger->info("delta_t: {}", delta_t);
    logger->info("Log Level: {}", log_level);

    logger->info("Force Type: {}",
                 (force_type == LennardJones ? "Lennard-Jones" : "Gravity"));

    if (container_type == LinkedCells) {
      using BoundaryType = LinkedCellsConfig::BoundaryType;
      logger->info("Container Type: Linked Cells");

      const auto& linked_cells_data =
          std::get<LinkedCellsConfig>(container_data);
      logger->info("-- Domain: ({}, {}, {})", linked_cells_data.domain[0],
                   linked_cells_data.domain[1], linked_cells_data.domain[2]);
      logger->info("-- Cutoff Radius: {}", linked_cells_data.cutoff_radius);

      logger->info("Boundary Configuration:");
      logger->info("------------------------");
      logger->info(
          "North Boundary: {}",
          (linked_cells_data.boundary_config.north == LinkedCellsConfig::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.north == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));
      logger->info(
          "South Boundary: {}",
          (linked_cells_data.boundary_config.south == BoundaryType::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.south == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));
      logger->info(
          "East Boundary: {}",
          (linked_cells_data.boundary_config.east == BoundaryType::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.east == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));
      logger->info(
          "West Boundary: {}",
          (linked_cells_data.boundary_config.west == BoundaryType::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.west == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));
      logger->info(
          "Up Boundary: {}",
          (linked_cells_data.boundary_config.up == BoundaryType::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.up == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));
      logger->info(
          "Down Boundary: {}",
          (linked_cells_data.boundary_config.down == BoundaryType::Outflow
               ? "Outflow"
           : linked_cells_data.boundary_config.down == BoundaryType::Reflective
               ? "Reflective"
               : "Periodic"));

    } else {
      logger->info("Container Type: Direct Sum");
    }

    logger->info("============================");
  }
};

#endif  // SIMULATION_H
