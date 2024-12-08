//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include <filesystem>

#include "debug/debug_print.h"
#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/types.h"
#include "forces/Gravity.h"
#include "forces/LennardJones.h"
#include "forces/SingularGravity.h"
#include "input.hxx"
#include "io/file/out/checkpoint-schema.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath,
                     Arguments& simulation_parameters) {
  const std::filesystem::path path(filepath);
  validate_path(path, ".xml", "Simulation input");
  try {
    const std::unique_ptr<::simulation> config = simulation_(filepath);
    SpdWrapper::get()->info("Reading XML file {}", filepath);
    auto& metadata = config->metadata();
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
      linked_cells_config.boundary_config = boundary_config;
      simulation_parameters.container_data = linked_cells_config;
      DEBUG_PRINT("Using LinkedCells container");
    } else {
      SpdWrapper::get()->warn(
          "No container provided, using default LinkedCells");
    }
    if (auto& force = metadata.force(); force.LennardJones().present()) {
      // auto _force = std::make_unique<LennardJones>();
      simulation_parameters.interactive_force_types.emplace_back(
          LennardJonesConfig{});
      DEBUG_PRINT("Using LennardJones");
    } else if (force.Gravity().present()) {
      // auto _force = std::make_unique<Gravity>();
      // simulation_parameters.interactive_forces.push_back(std::move(_force));
      simulation_parameters.interactive_force_types.emplace_back(
          GravityConfig{});
      DEBUG_PRINT("Using Gravity");
    } else {
      SpdWrapper::get()->warn("No force provided, using default LennardJones");
    }
    if (auto& singular_force = metadata.force();
        singular_force.SingularGravity().present()) {
      // auto _force = std::make_unique<SingularGravity>(
      //     singular_force.SingularGravity()->g().get());
      // simulation_parameters.singular_forces.push_back(std::move(_force));
      simulation_parameters.singular_force_types.emplace_back(
          SingularGravityConfig{singular_force.SingularGravity()->g().get()});
    }
    if (metadata.checkpoint().present()) {
      loadCheckpoint(metadata.checkpoint().get(), particles);
    }
    if (config->thermostat().present()) {
      auto thermostat = config->thermostat();
      ThermostatConfig thermostat_config = {
          .T_init = thermostat->T_init(),
          .T_target = thermostat->T_target(),  // TODO: make this optional
          .deltaT =
              std::numeric_limits<double>::max(),  // Default to infinity, dont
                                                   // use infinity() limit
                                                   // because of -ffast-math
          .n_thermostat = thermostat->n_thermostat(),
      };

      if (thermostat->deltaT().present()) {
        thermostat_config.deltaT = thermostat->deltaT().get();
      }

      SpdWrapper::get()->info("Checkpoint thermostat deltaT {}",
                              thermostat_config.deltaT);
      simulation_parameters.thermostat_config = thermostat_config;
      simulation_parameters.use_thermostat = true;
    } else {
      simulation_parameters.use_thermostat = false;
    }
    // TODO: singular forces here , but its a good start
    simulation_parameters.delta_t = metadata.delta_t();
    simulation_parameters.t_end = metadata.t_end();
    const auto& twoD = metadata.twoD();

    if (config->cuboids() != nullptr) {
      for (const auto& cubes : config->cuboids()->cuboid()) {
        const auto& _corner = cubes.corner();
        const auto& _dimensions = cubes.dimensions();
        const auto& _velocity = cubes.velocity();
        dvec3 corner = unwrapVec<const Dvec3Type&, dvec3>(_corner, "corner");
        ivec3 dimensions =
            unwrapVec<const Ivec3Type&, ivec3>(_dimensions, "dimensions");
        dvec3 velocity =
            unwrapVec<const Dvec3Type&, dvec3>(_velocity, "velocity");
        double mv;
        if (config->thermostat().present() && cubes.mv() == 0.0) {
          mv = std::sqrt(simulation_parameters.thermostat_config.T_init /
                         cubes.mass());
        } else {
          mv = cubes.mv();
        }

        CuboidGenerator cg(corner, dimensions, cubes.h(), cubes.mass(),
                           velocity, mv, cubes.epsilon(), cubes.sigma(),
                           cubes.type(), twoD);

        cg.generate(particles);
      }
    }

    if (config->spheroids() != nullptr) {
      for (const auto& spheres : config->spheroids()->spheroid()) {
        const auto& _origin = spheres.origin();
        const auto& _velocity = spheres.velocity();
        dvec3 origin = {_origin.x(), _origin.y(), _origin.z()};
        dvec3 velocity = {_velocity.x(), _velocity.y(), _velocity.z()};
        double mv;
        if (config->thermostat().present() && spheres.mv() == 0.0) {
          mv = std::sqrt(simulation_parameters.thermostat_config.T_init /
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
  SpdWrapper::get()->info("Looking and for checkpoint file");
  const std::filesystem::path filepath(_filepath);
  validate_path(filepath, ".xml", "checkpoint");
  try {
    DEBUG_PRINT("Found checkpoint file");
    const std::unique_ptr<::CheckpointType> checkpoint = Checkpoint(filepath);
    SpdWrapper::get()->info("Reading checkpoint particles");
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
      // TODO: this is also bad
      temp_particles.emplace_back(position, velocity, force, old_force, mass,
                                  epsilon, sigma, type);
    }
    SpdWrapper::get()->info("Read {} particles from checkpoint",
                            temp_particles.size());
    particles.reserve(particles.size() + temp_particles.size());
    particles.insert(particles.end(), temp_particles.begin(),
                     temp_particles.end());
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading checkpoint file: {}", e.what());
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