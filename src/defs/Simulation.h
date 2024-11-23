//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once
#include <stdexcept>

#include "defs/Particle.h"
#include "utils/SpdWrapper.h"

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
    logger->info("--------------------------");
    logger->info("Time End: {}", t_end);
    logger->info("Time Step: {}", delta_t);
    logger->info("Log Level: {}", log_level);

    logger->info("Force Type: {}",
                 (force_type == LennardJones ? "Lennard-Jones" : "Gravity"));

    logger->info("Container Type: {}",
                 (container_type == LinkedCells ? "Linked Cells" : "Direct Sum"));

    if (container_type == LinkedCells) {
      logger->info("Domain: ({}, {}, {})", domain[0], domain[1], domain[2]);
      logger->info("Cutoff Radius: {}", cutoff_radius);
    }

    logger->info("--------------------------");
  }
};

class Simulation {
public:

  void printConfiguration(Arguments& arguments);
};



// TODO: find a nice place for this
/**
 * @brief translates a vector from the xml parser to a valid "standard" c++ type
 * @tparam SVec source vector as one of the xsd types, e.g. IVec3Type
 * @tparam TVec target in c++-space, e.g ivec3 (:= std::array<int,3>)
 * @param source the vector passed form the xml parser
 * @param paramName the name of the parameter we are dealing with for throwing
 * an error
 * @return the vector in c++-space types.
 */
template <typename SVec, typename TVec>
TVec unwrapVec(const SVec& source, const std::string& paramName) {
  try {
    return TVec{source.x(), source.y(), source.z()};
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to unwrap vector " + paramName + ": " +
                             e.what());
  }
}

#endif  // SIMULATION_H
