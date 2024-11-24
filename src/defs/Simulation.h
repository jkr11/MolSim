//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include "defs/types.h"
#include "utils/SpdWrapper.h"

/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  double t_end;
  double delta_t;
  // double output_time_step_size; // TODO: this is still needed, this is
  // supposed to be fps later on
  std::string log_level;
  enum ForceType { LennardJones, Gravity } force_type;
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
  enum ContainerType { LinkedCells, DirectSum } container_type;
  // TODO: make this into an std::variant<LinkedCellsStruct,DirectSumStruct> as
  // we add more information that directsum doesnt need to know about.
  void printConfiguration() const {
    const auto logger = SpdWrapper::get();

    logger->info("Simulation Configuration:");
    logger->info("============================");
    logger->info("t_end: {}", t_end);
    logger->info("delta_t: {}", delta_t);
    logger->info("Log Level: {}", log_level);

    logger->info("Force Type: {}",
                 (force_type == LennardJones ? "Lennard-Jones" : "Gravity"));

    logger->info(
        "Container Type: {}",
        (container_type == LinkedCells ? "Linked Cells" : "Direct Sum"));

    if (container_type == LinkedCells) {
      logger->info("-- Domain: ({}, {}, {})", domain[0], domain[1], domain[2]);
      logger->info("-- Cutoff Radius: {}", cutoff_radius);
      logger->info("-- Boundary Configuration:");
      logger->info("---- North Boundary: {}",
                   (boundary_config.north == Outflow ? "Outflow" :
                    boundary_config.north == Reflective ? "Reflective" : "Periodic"));
      logger->info("---- South Boundary: {}",
                   (boundary_config.south == Outflow ? "Outflow" :
                    boundary_config.south == Reflective ? "Reflective" : "Periodic"));
      logger->info("---- East Boundary: {}",
                   (boundary_config.east == Outflow ? "Outflow" :
                    boundary_config.east == Reflective ? "Reflective" : "Periodic"));
      logger->info("---- West Boundary: {}",
                   (boundary_config.west == Outflow ? "Outflow" :
                    boundary_config.west == Reflective ? "Reflective" : "Periodic"));
      logger->info("---- Up Boundary: {}",
                   (boundary_config.up == Outflow ? "Outflow" :
                    boundary_config.up == Reflective ? "Reflective" : "Periodic"));
      logger->info("---- Down Boundary: {}",
                   (boundary_config.down == Outflow ? "Outflow" :
                    boundary_config.down == Reflective ? "Reflective" : "Periodic"));
    }

    logger->info("============================");
  }
};

#endif  // SIMULATION_H
