//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "defs/types.h"
#include "input.hxx"
#include "spdlog/fmt/bundled/args.h"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath) {
  try {
    const std::unique_ptr<::simulation> config = simulation_(filepath);
    SpdWrapper::get()->info("Reading XML file {}", filepath);
    if (config->metadata() != nullptr) {
      auto& metadata = config->metadata();
      if (metadata->delta_t() != nullptr) {
        simulation_parameters.delta_t = metadata->delta_t().get();
      }
      if (metadata->t_end() != nullptr) {
        simulation_parameters.t_end = metadata->t_end().get();
      }
      if (metadata->r_cutoff() != nullptr) {
        simulation_parameters.cutoff_radius = metadata->r_cutoff().get();
      }
      if (metadata->domain() != nullptr) {
        simulation_parameters.domain =
            unwrapVec<Ivec3Type&, ivec3>(metadata->domain().get(), "domain");
      }
    }

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
