//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include <filesystem>
#include <iostream>

#include "debug/debug_print.h"
#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/MembraneGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/types.h"
#include "forces/Gravity.h"
#include "input.hxx"
#include "io/file/out/checkpoint-schema.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath,
                     Arguments& simulation_parameters) {
  const std::filesystem::path path(filepath);
  validatePath(path, ".xml", "Simulation input");
  try {
    const std::unique_ptr<::simulation> config = simulation_(filepath);
    INFO_FMT("Reading XML file {}", filepath);
    auto& metadata = config->metadata();
    INFO("Reading containers ...")
    if (auto& container = metadata.container();
        container.directSum().present()) {
      DirectSumConfig direct_sum_config;
      simulation_parameters.container_data = direct_sum_config;
      DEBUG_PRINT("Using DirectSum container")
    } else if (container.linkedCells().present()) {
      LinkedCellsConfig linked_cells_config{};
      const auto& domain = container.linkedCells().get().domain();
      linked_cells_config.domain =
          unwrapVec<const Ivec3Type&, ivec3>(domain, "domain");
      linked_cells_config.cutoff_radius =
          container.linkedCells().get().r_cutoff();
      const auto& boundaries = container.linkedCells().get().boundary();
      const LinkedCellsConfig::BoundaryConfig boundary_config = {
          .x_high = toBoundaryType(boundaries.x_high()),
          .x_low = toBoundaryType(boundaries.x_low()),
          .y_high = toBoundaryType(boundaries.y_high()),
          .y_low = toBoundaryType(boundaries.y_low()),
          .z_high = toBoundaryType(boundaries.z_high()),
          .z_low = toBoundaryType(boundaries.z_low()),
      };
      validateBoundaries(boundary_config);

      linked_cells_config.boundary_config = boundary_config;
      simulation_parameters.container_data = linked_cells_config;
      DEBUG_PRINT("Using LinkedCells container");
    } else {
      SpdWrapper::get()->warn(
          "No container provided, using default LinkedCells");
    }
    INFO("Reading interactive forces ...")
    if (auto& force = metadata.force(); force.LennardJones().present()) {
      simulation_parameters.interactive_force_types.emplace_back(
          LennardJonesConfig{});
      DEBUG_PRINT("Using LennardJones");
    } else if (force.Gravity().present()) {
      simulation_parameters.interactive_force_types.emplace_back(
          GravityConfig{});
      DEBUG_PRINT("Using Gravity");
    } else if (force.TruncatedLennardJonesForce().present()) {
      simulation_parameters.interactive_force_types.emplace_back(
          TruncatedLennardJonesConfig{});
      DEBUG_PRINT("Using TruncatedLennardJonesForce");
    } else {
      SpdWrapper::get()->warn("No force provided, using default LennardJones");
    }
    INFO("Reading singular foces ...")
    if (auto& singular_force = metadata.force();
        singular_force.SingularGravity().present()) {
      SpdWrapper::get()->info("Adding singular gravity on axis {}",
                              singular_force.SingularGravity()->axis());

      simulation_parameters.singular_force_types.emplace_back(
          SingularGravityConfig{singular_force.SingularGravity()->g(),
                                singular_force.SingularGravity()->axis()});
    }
    if (auto& singular_force = metadata.force();
        singular_force.HarmonicForce().present()) {
      simulation_parameters.singular_force_types.emplace_back(
          HarmonicForceConfig{singular_force.HarmonicForce()->r_0(),
                              singular_force.HarmonicForce()->k()});
      DEBUG_PRINT("Using HarmonicForce");
    }
    if (auto& index_force = metadata.force();
        index_force.IndexForce().present()) {
      SpdWrapper::get()->info("Building index force");
      std::vector<ivec3> indeces{};
      for (auto& i : index_force.IndexForce()->index()) {
        indeces.push_back(unwrapVec<const Ivec3Type&, ivec3>(i, "index"));
      }
      SpdWrapper::get()->info("indeces[0] {}", indeces[0][0]);
      IndexForceConfig index_force_config{
          indeces,
          {},
          index_force.IndexForce()->time(),
          unwrapVec<const Dvec3Type&, dvec3>(
              index_force.IndexForce()->force_values(), "force_values"),
      };
      simulation_parameters.index_force_configs.push_back(index_force_config);
    }
    INFO("Checking if checkpoint is  present ...")
    if (metadata.checkpoint().present()) {
      INFO("Membrane checkpoint present")
      auto name = static_cast<std::string>(metadata.checkpoint()->name());
      if (bool is_membrane = metadata.checkpoint()->is_membrane()) {
        ivec3 dimensions = unwrapVec<Ivec3Type, ivec3>(
            metadata.checkpoint()->domain(), "domain of membrane checkpoint");
        loadCheckpointMembrane(name, dimensions, particles);
      } else {
        loadCheckpoint(name, particles);
      }
    }

    StatisticsConfig statistics_config;
    if (metadata.statistics().present()) {
      auto statistics = metadata.statistics().get();
      statistics_config = {
          .calc_stats = true,
          .x_bins = statistics.x_bins(),
          .y_bins = statistics.y_bins(),
          .output_interval = statistics.output_interval(),
          .velocity_output_location = "velocity.csv",
          .density_output_location = "density.csv",
      };
    } else {
      statistics_config = {
          .calc_stats = false,
          .velocity_output_location = "velocity.csv",
          .density_output_location = "density.csv",
      };
    }
    if (metadata.use_c18_strategy().present()) {
      if (metadata.use_c18_strategy().get()) {
        INFO("Using c18 strategy")
        simulation_parameters.strategy = ParallelStrategy::STRATEGY_2;
      } else {
        INFO("Using force buffers reader")
        simulation_parameters.strategy = ParallelStrategy::STRATEGY_1;
      }
    } else {
      simulation_parameters.strategy = ParallelStrategy::STRATEGY_3;
    }
    simulation_parameters.statistics_config = statistics_config;

    INFO("Checking for thermostat ...")
    if (config->thermostat().present()) {
      auto thermostat = config->thermostat();
      ThermostatConfig thermostat_config = {
          .t_init = thermostat->T_init(),
          .t_target =
              thermostat->T_init(),  // we initialize both this and deltaT to
                                     // its defaults first and then check if the
                                     // actual values are present
          .delta_t =
              std::numeric_limits<double>::max(),  // Default to infinity, dont
                                                   // use infinity() limit
                                                   // because of -ffast-math
          .n_thermostat = thermostat->n_thermostat(),
          .use_thermal_motion =
              static_cast<bool>(thermostat->use_thermal_motion()),
          .two_d = static_cast<bool>(config->metadata().twoD()),
      };

      if (thermostat->deltaT().present()) {
        thermostat_config.delta_t = thermostat->deltaT().get();
      }
      if (thermostat->T_target().present()) {
        thermostat_config.t_target = thermostat->T_target().get();
      }

      simulation_parameters.thermostat_config = thermostat_config;
      simulation_parameters.use_thermostat = true;
    } else {
      simulation_parameters.use_thermostat = false;
    }
    simulation_parameters.delta_t = metadata.delta_t();
    simulation_parameters.t_end = metadata.t_end();

    const auto& twoD = metadata.twoD();
    INFO_FMT("Simulation has dimension {}", twoD ? 2 : 3)
    INFO("Reading in particle shape configurations ...")
    if (config->cuboids() != nullptr) {
      for (const auto& cubes : config->cuboids()->cuboid()) {
        SpdWrapper::get()->info("Generating cuboid");
        const auto& _corner = cubes.corner();
        const auto& _dimensions = cubes.dimensions();
        const auto& _velocity = cubes.velocity();
        dvec3 corner = unwrapVec<const Dvec3Type&, dvec3>(_corner, "corner");
        ivec3 dimensions =
            unwrapVec<const Ivec3Type&, ivec3>(_dimensions, "dimensions");
        dvec3 velocity =
            unwrapVec<const Dvec3Type&, dvec3>(_velocity, "velocity");
        double mv;
        if (config->thermostat().present() &&
            velocity == dvec3{0.0, 0.0, 0.0}) {
          mv = std::sqrt(simulation_parameters.thermostat_config.t_init /
                         cubes.mass());
        } else {
          mv = cubes.mv();
        }

        // These info print themselves
        CuboidGenerator cg(corner, dimensions, cubes.h(), cubes.mass(),
                           velocity, mv, cubes.epsilon(), cubes.sigma(),
                           cubes.type(), twoD);
        cg.generate(particles);
      }
    }

    if (config->membranes() != nullptr) {
      std::cout << "membrane_config found" << std::endl;
      for (const auto& membranes : config->membranes()->membrane()) {
        SpdWrapper::get()->info("Generating membranes");
        const auto& _corner = membranes.corner();
        const auto& _dimensions = membranes.dimensions();
        const auto& _velocity = membranes.velocity();
        dvec3 corner = unwrapVec<const Dvec3Type&, dvec3>(_corner, "corner");
        ivec3 dimensions =
            unwrapVec<const Ivec3Type&, ivec3>(_dimensions, "dimensions");
        dvec3 velocity =
            unwrapVec<const Dvec3Type&, dvec3>(_velocity, "velocity");
        double mv;
        if (config->thermostat().present() &&
            velocity == dvec3{0.0, 0.0, 0.0}) {
          mv = std::sqrt(simulation_parameters.thermostat_config.t_init /
                         membranes.mass());
        } else {
          mv = membranes.mv();
        }
        MembraneGenerator mg(
            corner, dimensions, membranes.h(), membranes.mass(), velocity, mv,
            membranes.epsilon(), membranes.sigma(), membranes.type(), twoD,
            simulation_parameters.index_force_configs[0].indeces);
        mg.generate(particles);
        std::vector<int> ids = mg.getIndices();
        simulation_parameters.index_force_configs[0].ids = ids;
        if (std::holds_alternative<LinkedCellsConfig>(
                simulation_parameters.container_data)) {
          LinkedCellsConfig linked_cells_config_tmp =
              std::get<LinkedCellsConfig>(simulation_parameters.container_data);
          linked_cells_config_tmp.is_membrane = true;
          simulation_parameters.container_data = linked_cells_config_tmp;
        }
      }
    }

    if (config->spheroids() != nullptr) {
      for (const auto& spheres : config->spheroids()->spheroid()) {
        const auto& _origin = spheres.origin();
        const auto& _velocity = spheres.velocity();
        dvec3 origin = {_origin.x(), _origin.y(), _origin.z()};
        dvec3 velocity = {_velocity.x(), _velocity.y(), _velocity.z()};
        double mv;
        if (config->thermostat().present() &&
            velocity == dvec3{0.0, 0.0, 0.0}) {
          mv = std::sqrt(simulation_parameters.thermostat_config.t_init /
                         spheres.mass());
        } else {
          mv = spheres.mv();
        }

        SpheroidGenerator sg(origin, spheres.radius(), spheres.h(),
                             spheres.mass(), velocity, spheres.epsilon(),
                             spheres.sigma(), spheres.type(), mv, twoD);

        sg.generate(particles);
      }
    }
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
    exit(EXIT_FAILURE);
  }
}

void XmlReader::loadCheckpoint(const std::string& _filepath,
                               std::vector<Particle>& particles) {
  const std::filesystem::path filepath(_filepath);
  validatePath(filepath, ".xml", "checkpoint");
  try {
    DEBUG_PRINT("Found checkpoint file");
    const std::unique_ptr<::CheckpointType> checkpoint = Checkpoint(filepath);
    INFO("Reading checkpoint particles");
    std::vector<Particle> temp_particles;
    for (const auto& p : checkpoint->Particles().Particle()) {
      auto position =
          unwrapVec<const CDvec3Type, dvec3>(p.Position(), "position");
      auto velocity =
          unwrapVec<const CDvec3Type, dvec3>(p.Velocity(), "position");
      auto force = unwrapVec<const CDvec3Type, dvec3>(p.Force(), "force");
      auto old_force =
          unwrapVec<const CDvec3Type, dvec3>(p.OldForce(), "old_force");
      double mass = p.mass();
      double epsilon = p.epsilon();
      double sigma = p.sigma();
      int type = p.type();
      temp_particles.emplace_back(position, velocity, force, old_force, mass,
                                  type, epsilon, sigma);
    }
    INFO_FMT("Read {} particles from checkpoint", temp_particles.size());
    particles.reserve(particles.size() + temp_particles.size());
    particles.insert(particles.end(), temp_particles.begin(),
                     temp_particles.end());
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading checkpoint file: {}", e.what());
  }
}

void XmlReader::loadCheckpointMembrane(const std::string& _filepath,
                                       const ivec3& dimensions,
                                       std::vector<Particle>& particles) {
  INFO("Loading membrane checkpoint")
  const std::filesystem::path filepath(_filepath);
  validatePath(filepath, ".xml", "checkpoint");
  try {
    DEBUG_PRINT("Found checkpoint file");
    const std::unique_ptr<::CheckpointType> checkpoint = Checkpoint(filepath);
    INFO("Reading checkpoint particles");
    std::vector<Particle> temp_particles;
    for (const auto& p : checkpoint->Particles().Particle()) {
      auto position =
          unwrapVec<const CDvec3Type, dvec3>(p.Position(), "position");
      auto velocity =
          unwrapVec<const CDvec3Type, dvec3>(p.Velocity(), "position");
      auto force = unwrapVec<const CDvec3Type, dvec3>(p.Force(), "force");
      auto old_force =
          unwrapVec<const CDvec3Type, dvec3>(p.OldForce(), "old_force");
      double mass = p.mass();
      double epsilon = p.epsilon();
      double sigma = p.sigma();
      int type = p.type();
      temp_particles.emplace_back(position, velocity, force, old_force, mass,
                                  type, epsilon, sigma);
    }
    INFO_FMT("Read {} particles from checkpoint", temp_particles.size());
    std::size_t previous_size = particles.size();
    particles.reserve(previous_size + temp_particles.size());
    particles.insert(particles.end(), temp_particles.begin(),
                     temp_particles.end());

    for (int i = 0; i < dimensions[0]; i++) {
      for (int j = 0; j < dimensions[1]; j++) {
        for (int k = 0; k < dimensions[2]; k++) {
          const std::size_t current_index =
              i * dimensions[1] * dimensions[2] + j * dimensions[2] + k;

          // Iterate over all neighbors including diagonals
          for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
              for (int dk = -1; dk <= 1; dk++) {
                if (di == 0 && dj == 0 && dk == 0) {
                  continue;
                }

                const long ni = i + di;
                const long nj = j + dj;
                const long nk = k + dk;

                if (ni >= 0 && ni < dimensions[0] && nj >= 0 &&
                    nj < dimensions[1] && nk >= 0 && nk < dimensions[2]) {
                  const long neighbor_index =
                      ni * dimensions[1] * dimensions[2] + nj * dimensions[2] +
                      nk;
                  const bool is_diagonal =
                      (di != 0) + (dj != 0) + (dk != 0) > 1;
                  Particle* neighbor_particle =
                      &particles[previous_size + neighbor_index];
                  particles[previous_size + current_index].pushBackNeighbour(
                      is_diagonal, reinterpret_cast<size_t>(neighbor_particle));
                }
              }
            }
          }
        }
      }
    }
    DEBUG_PRINT("particles: " + std::to_string(particles.size()));

  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading checkpoint file: {}", e.what());
  }
}

void XmlReader::validateBoundaries(
    const LinkedCellsConfig::BoundaryConfig& boundary) {
  if (boundary.x_high != boundary.x_low &&
      (boundary.x_high == LinkedCellsConfig::Periodic ||
       boundary.x_low == LinkedCellsConfig::Periodic)) {
    throw std::runtime_error("x dimension has incompatible boundaries");
  }

  if (boundary.y_high != boundary.y_low &&
      (boundary.y_high == LinkedCellsConfig::Periodic ||
       boundary.y_low == LinkedCellsConfig::Periodic)) {
    throw std::runtime_error("y dimension has incompatible boundaries");
  }

  if (boundary.z_high != boundary.z_low &&
      (boundary.z_high == LinkedCellsConfig::Periodic ||
       boundary.z_low == LinkedCellsConfig::Periodic)) {
    throw std::runtime_error("z dimension has incompatible boundaries");
  }
}

using LBoundaryType =
    LinkedCellsConfig::BoundaryType;  // BoundaryType is double used by the
                                      // xsd config files so be careful
template <typename BT>
LBoundaryType toBoundaryType(const BT& boundary_type) {
  if (boundary_type.Outflow().present()) {
    return LBoundaryType::Outflow;
  }
  if (boundary_type.Periodic().present()) {
    return LBoundaryType::Periodic;
  }
  if (boundary_type.Reflective().present()) {
    return LBoundaryType::Reflective;
  }
  throw std::runtime_error("Unknown boundary type");
}

template <typename T>
double unpackOpt(const T& opt) {
  if (opt.present()) {
    return opt.value();
  }
  return std::numeric_limits<double>::infinity();
}