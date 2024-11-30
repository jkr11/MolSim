//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <variant>

#include "Generators/ParticleGenerator.h"
#include "calc/VerletIntegrator.h"
#include "defs/types.h"
#include "forces/Force.h"
#include "utils/SpdWrapper.h"

/**
 * @brief holds the specification for the LinkedCellsContainer
 */
struct LinkedCellsConfig {
  ivec3 domain;
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
 * @brief input for the CuboidGenerator class
 */
struct CuboidConfig {
  dvec3 left_corner;
  ivec3 dimensions;
  dvec3 initial_velocity;
  double h;
  double mass;
  double mean_velocity;
  double epsilon;
  double sigma;
  int type;
  bool twoD;
};

/**
 * @brief input for the sphereoidGenerator class
 */
struct SphereoidConfig {
  dvec3 origin;
  double radius;
  dvec3 initial_velocity;
  double h;
  double mass;
  double mean_velocity;
  double epsilon;
  double sigma;
  int type;
  bool twoD;
};

/**
 * @brief extensible type for holding the input shapes
 */
using GeneratorConfig = std::variant<CuboidConfig, SphereoidConfig>;

/**
 * @brief struct to hold command line arguments
 */
struct Arguments {
  double t_end;
  double delta_t;
  enum ForceType { LennardJones, Gravity } force_type;
  std::variant<LinkedCellsConfig, DirectSumConfig> container_data;
  std::vector<std::unique_ptr<ParticleGenerator>> generator_configs;
};

/**
 * This class holds all the necessary info and objects to run the simulation.
 * This is done to clean up main and choosing a place to allocate the particle
 * vector.
 */
class Simulation {
 private:
  Arguments arguments;
  std::vector<Particle> particles;
  std::unique_ptr<ParticleContainer> container;
  std::unique_ptr<Force> forces;
  std::unique_ptr<VerletIntegrator> integrator;
  std::string output_directory;
  double step_size;

 public:
  explicit Simulation(Arguments arguments, std::string output_directory,
                      double step_size);
  ~Simulation();

  /**
   * initializes all the necessary objects such as the containers and forces.
   */
  void initParticles();
  /**
   * unpacks the specification from the generators into the particle vector.
   */
  void initParams();
  /**
   * runs the main loop of the program.
   */
  void run() const;
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

    const auto& linked_cells_data =
        std::get<LinkedCellsConfig>(args.container_data);
    logger->info("-- Domain: ({}, {}, {})", linked_cells_data.domain[0],
                 linked_cells_data.domain[1], linked_cells_data.domain[2]);
    logger->info("-- Cutoff Radius: {}", linked_cells_data.cutoff_radius);

    logger->info("Boundary Configuration:");
    logger->info("------------------------");

    const auto& [x_high, x_low, y_high, y_low, z_high, z_low] =
        linked_cells_data.boundary_config;
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
