//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once

#include <variant>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "defs/types.h"
#include "forces/SingularForce.h"
#include "utils/SpdWrapper.h"

/**
 * @brief Parallelization Strategy
 */
enum ParallelStrategy { STRATEGY_1, STRATEGY_2, STRATEGY_3 };

/**
 * @brief Index Force config struct
 */
struct IndexForceConfig {
  /**
   * Cuboid coordinates to be targetted
   */
  std::vector<ivec3> indeces{};
  /**
   * Indices to be targetted
   */
  std::vector<int> ids{};

  /**
   * End time of the force application
   */
  double end_time{};
  /**
   * force vector applied
   */
  dvec3 force_values{};
};

/**
 * @brief holds the specification for the LinkedCellsContainer
 */
struct LinkedCellsConfig {
  /**
   * Domain size
   */
  ivec3 domain;
  /**
   * Cutoff radius of the Linked-Cells algorithm
   */
  double cutoff_radius;
  /**
   * Boundary Type enum for different boundary conditions
   */
  enum BoundaryType { Outflow, Reflective, Periodic } boundary_type;
  /**
   * Struct which holds the boundary config for the linked cells container
   */
  struct BoundaryConfig {
    /**
     * high x boundary type
     */
    BoundaryType x_high;
    /**
     * low x boundary type
     */
    BoundaryType x_low;
    /**
     * high y boundary type
     */
    BoundaryType y_high;
    /**
     * low y boundary type
     */
    BoundaryType y_low;
    /**
     * high z boundary type
     */
    BoundaryType z_high;
    /**
     * low z boundary type
     */
    BoundaryType z_low;
  } boundary_config;

  /**
   * Whether the simulation uses membranes or not
   */
  bool is_membrane = false;
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
  /**
   * "gravitational" force applied
   */
  double g{};
  /**
   * axis in which the force is applied
   */
  int axis{};
};

/**
 * @brief holds the harmonic force parameters
 */
struct HarmonicForceConfig {
  /**
   * average bond length
   */
  double r_0{};
  /**
   * Spring/stiffness constant
   */
  double k{};
};

/**
 * @brief holds the LennardJones force parameters
 */
struct LennardJonesConfig {};

/**
 * @brief hold the Truncated LennardJones force parameters
 */
struct TruncatedLennardJonesConfig {};

/**
 * @brief holds the interactive gravity parameters
 */
struct GravityConfig {};

/**
 * @brief holds instance data for Thermostat
 */
struct ThermostatConfig {
  /**
   * init temperature
   */
  double t_init{};
  /**
   * target temperature
   */
  double t_target{};
  /**
   * delta time
   */
  double delta_t{};
  /**
   * periodicity of application
   */
  int n_thermostat{};
  /**
   * Whether thermostat is relative
   */
  bool use_relative{};
  /**
   * Whether thermostat uses thermal motion
   */
  bool use_thermal_motion{};
  /**
   * Whether thermostat is 2D or 3D
   */
  bool two_d{};
};

/**
 * @brief Statistics configuration struct
 */
struct StatisticsConfig {
  /**
   * TODO ?
   */
  bool calc_stats{};
  /**
   * Number of x bins
   */
  int x_bins{};
  /**
   * Number of y bins
   */
  int y_bins{};
  /**
   * TODO ?
   */
  int output_interval{};
  /**
   * TODO ?
   */
  std::string velocity_output_location{};
  /**
   * TODO ?
   */
  std::string density_output_location{};
};

//TODO: delete?
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

//TODO: delete?
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

//TODO: delete?
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

/**
 * @brief Struct which hold the simulation arguments
 */
struct Arguments {
  /**
   * simulation end time
   */
  double t_end;
  /**
   * simulation delta time
   */
  double delta_t;
  /**
   * Thermostat config
   */
  ThermostatConfig thermostat_config;
  /**
   * Whether to use the thermostat
   */
  bool use_thermostat;
  /**
   * Which container to use
   */
  std::variant<LinkedCellsConfig, DirectSumConfig> container_data;
  /**
   * Vector of applied singular forces
   */
  std::vector<SingularForceTypes> singular_force_types;
  /**
   * Vector of applied interactive forces
   */
  std::vector<InteractiveForceTypes> interactive_force_types;
  /**
   * Vector of applied index forces
   */
  std::vector<IndexForceConfig> index_force_configs;
  /**
   * Statistics configuration
   */
  StatisticsConfig statistics_config;

  //TODO
  SphereoidGeneratorConfig spheroid_generator_config;
  //TODO
  CuboidGeneratorConfig cuboid_generator_config;
  //TODO
  MembraneGeneratorConfig membrane_generator_config;

  /**
   * Parallelization strategy
   */
  ParallelStrategy strategy;
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

    const auto& [domain, cutoff_radius, boundary_type, boundary_config,
                 is_membrane] =
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
#ifdef _OPENMP
  logger->info("Number of Threads: {}", omp_get_max_threads());
#endif
  logger->info("============================");
}

#endif  // SIMULATION_H
