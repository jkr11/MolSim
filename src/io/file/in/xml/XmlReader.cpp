//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/Simulation.h"
#include "input.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath) {
  try {
    const std::unique_ptr<::simulation> config = simulation_(filepath);
    SpdWrapper::get()->info("Reading XML file {}", filepath);
    auto& metadata = config->metadata();
    if (auto& container = metadata.container();
        container.directSum().present()) {
      simulation_parameters.container_type = Arguments::DirectSum;
    } else if (container.linkedCells().present()) {
      simulation_parameters.container_type = Arguments::LinkedCells;
      const auto& domain = container.linkedCells().get().domain();
      simulation_parameters.domain = unwrapVec<const Ivec3Type&, ivec3>(domain, "domain");
      simulation_parameters.cutoff_radius = container.linkedCells().get().r_cutoff();
    } else {
      SpdWrapper::get()->info(
          "No container provided, using default LinkedCells");
    }
    simulation_parameters.delta_t = metadata.delta_t();
    simulation_parameters.t_end = metadata.t_end();

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
                           cubes.type(), cubes.twoD());

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
                             spheres.sigma(), spheres.type(), spheres.twoD());

        sg.generate(particles);
      }
    }

  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
  }
}

std::tuple<double, double, double, ivec3> XmlReader::pass() const {
  return {simulation_parameters.delta_t, simulation_parameters.t_end,
          simulation_parameters.cutoff_radius, simulation_parameters.domain};
}
