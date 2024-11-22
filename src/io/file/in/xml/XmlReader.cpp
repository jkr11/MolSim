//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include <filesystem>

#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/Simulation.h"
#include "input.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath) {
  const std::filesystem::path path(filepath);
  if (!exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }
  if (path.extension() != ".xml") {
    throw std::invalid_argument("File extension is not supported: " +
                                path.string());
  }
  try {
    const std::unique_ptr<::simulation> config = simulation_(filepath);
    SpdWrapper::get()->info("Reading XML file {}", filepath);
    auto& metadata = config->metadata();
    if (auto& container = metadata.container();
        container.directSum().present()) {
      simulation_parameters.container_type = Arguments::DirectSum;
      SpdWrapper::get()->info("Using directSumContainer");
    } else if (container.linkedCells().present()) {
      simulation_parameters.container_type = Arguments::LinkedCells;
      const auto& domain = container.linkedCells().get().domain();
      simulation_parameters.domain =
          unwrapVec<const Ivec3Type&, ivec3>(domain, "domain");
      simulation_parameters.cutoff_radius =
          container.linkedCells().get().r_cutoff();
      SpdWrapper::get()->info("Using LinkedCellsContainer");
    } else {
      SpdWrapper::get()->info(
          "No container provided, using default LinkedCells");
    }
    if (auto& force = metadata.force(); force.LennardJones().present()) {
      simulation_parameters.force_type = Arguments::LennardJones;
    } else if (force.Gravity().present()) {
      simulation_parameters.force_type = Arguments::Gravity;
    }
    simulation_parameters.delta_t = metadata.delta_t();
    simulation_parameters.t_end = metadata.t_end();
    // it's enough for twoD to be a global parameter
    // I didn't really see any need for a single simulation to support both twoD
    // and threeD interactions
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

        CuboidGenerator cg(corner, dimensions, cubes.h(), cubes.mass(),
                           velocity, cubes.mv(), cubes.epsilon(), cubes.sigma(),
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

        SpheroidGenerator sg(origin, spheres.radius(), spheres.h(),
                             spheres.mass(), velocity, spheres.epsilon(),
                             spheres.sigma(), spheres.type(), twoD);

        sg.generate(particles);
      }
    }

  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
  }
}

std::tuple<double, double, double, ivec3, Arguments::ForceType,
           Arguments::ContainerType>
XmlReader::pass() const {
  return {simulation_parameters.delta_t,
          simulation_parameters.t_end,
          simulation_parameters.cutoff_radius,
          simulation_parameters.domain,
          simulation_parameters.force_type,
          simulation_parameters.container_type};
}
