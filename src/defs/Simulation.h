//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once
#include <stdexcept>

#include "utils/SpdWrapper.h"

using dvec3 = std::array<double, 3>;
using ivec3 = std::array<int, 3>;

/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  // std::string input_file;
  double t_end;
  double delta_t;
  // double output_time_step_size; // TODO: this is still needed, this is
  // supposed to be fps later on
  std::string log_level;
  enum ForceType { LennardJones, Gravity } force_type;
  ivec3 domain;
  double cutoff_radius;
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

    logger->info("Container Type: {}",
                 (container_type == LinkedCells ? "Linked Cells" : "Direct Sum"));

    if (container_type == LinkedCells) {
      logger->info("-- Domain: ({}, {}, {})", domain[0], domain[1], domain[2]);
      logger->info("-- Cutoff Radius: {}", cutoff_radius);
    }

    logger->info("============================");
  }
};

// TODO: find a nice place for this


#endif  // SIMULATION_H
