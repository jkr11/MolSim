//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include <filesystem>

#include "debug/debug_print.h"
#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/types.h"
#include "input.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(const std::string& filepath,
                     Arguments& simulation_parameters) {
  const std::filesystem::path path(filepath);
  if (!exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }
  if (path.extension() != ".xml") {
    throw std::invalid_argument("File extension is not supported: " +
                                path.string());
  }
  try {
    const std::unique_ptr<simulation> config = simulation_(filepath);
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
      simulation_parameters.force_type = Arguments::LennardJones;
      DEBUG_PRINT("Using LennardJones");
    } else if (force.Gravity().present()) {
      simulation_parameters.force_type = Arguments::Gravity;
      DEBUG_PRINT("Using Gravity");
    } else {
      SpdWrapper::get()->warn("No force provided, using default LennardJones");
    }
    simulation_parameters.delta_t = metadata.delta_t();
    simulation_parameters.t_end = metadata.t_end();
    const auto& twoD = metadata.twoD();

    if (config->cuboids() != nullptr) {
      for (const auto& cubes : config->cuboids()->cuboid()) {
        const auto& _corner = cubes.corner();
        const auto& _dimensions = cubes.dimensions();
        const auto& _velocity = cubes.velocity();
        const dvec3 corner =
            unwrapVec<const Dvec3Type&, dvec3>(_corner, "corner");
        const ivec3 dimensions =
            unwrapVec<const Ivec3Type&, ivec3>(_dimensions, "dimensions");
        const dvec3 velocity =
            unwrapVec<const Dvec3Type&, dvec3>(_velocity, "velocity");
        /*
        CuboidConfig cuboid_config = {.left_corner = corner,
                                      .dimensions = dimensions,
                                      .initial_velocity = velocity,
                                      .h = cubes.h(),
                                      .mass = cubes.mass(),
                                      .mean_velocity = cubes.mv(),
                                      .epsilon = cubes.epsilon(),
                                      .sigma = cubes.sigma(),
                                      .type = cubes.type(),
                                      .twoD = twoD};
        simulation_parameters.generator_configs.emplace_back(cuboid_config);
        */
        CuboidGenerator cg(corner, dimensions, cubes.h(), cubes.mass(),
                           velocity, cubes.mv(), cubes.epsilon(), cubes.sigma(),
                           cubes.type(), twoD);
        simulation_parameters.generator_configs.emplace_back(
            std::make_unique<CuboidGenerator>(std::move(cg)));

        //simulation_parameters.generator_configs.emplace_back(std::make_unique<CuboidGenerator>(corner, dimensions, cubes.h(), cubes.mass(),
        //           velocity, cubes.mv(), cubes.epsilon(), cubes.sigma(),
        //           cubes.type(), twoD));

        // cg.generate(particles);
      }
    }

    if (config->spheroids() != nullptr) {
      for (const auto& spheres : config->spheroids()->spheroid()) {
        const auto& _origin = spheres.origin();
        const auto& _velocity = spheres.velocity();
        const dvec3 origin = {_origin.x(), _origin.y(), _origin.z()};
        const dvec3 velocity = {_velocity.x(), _velocity.y(), _velocity.z()};

        /*
        SphereoidConfig sphereoid_config = {
            .origin = origin,
            .radius = static_cast<double>(spheres.radius()),
            .initial_velocity = velocity,
            .h = spheres.h(),
            .mass = spheres.mass(),
            .epsilon = spheres.epsilon(),
            .sigma = spheres.sigma(),
            .type = spheres.type(),
            .twoD = twoD};

        simulation_parameters.generator_configs.emplace_back(spheroid_config);
        */
        SpheroidGenerator sg(origin, spheres.radius(), spheres.h(),
                             spheres.mass(), velocity, spheres.epsilon(),
                             spheres.sigma(), spheres.type(), twoD);
        simulation_parameters.generator_configs.emplace_back(
            std::make_unique<SpheroidGenerator>(std::move(sg)));

        // sg.generate(particles);
      }
    }
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
    exit(EXIT_FAILURE);
  }
}

using LBoundaryType =
    LinkedCellsConfig::BoundaryType;  // BoundaryType is double used by the xsd
                                      // config files so be careful
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