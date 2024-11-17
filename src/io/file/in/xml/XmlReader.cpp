//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include "defs/Generators/CuboidGenerator.h"
#include "defs/Generators/SpheroidGenerator.h"
#include "input.hxx"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath) {
  try {
    const std::unique_ptr< ::simulation> config = simulation_(filepath);
    SpdWrapper::get()->info("Reading XML file {}", filepath);
    for (const auto& cubes : config->cuboids().cuboid()) {
      const auto& _corner = cubes.corner();
      const auto& _dimensions = cubes.dimensions();
      const auto& _velocity = cubes.velocity();
      dvec3 corner = {_corner.x(), _corner.y(), _corner.z()};
      std::array<int, 3> dimensions = {_dimensions.x(), _dimensions.y(),
                                       _dimensions.z()};
      dvec3 velocity = {_velocity.x(), _velocity.y(), _velocity.z()};
      CuboidGenerator cg(corner, dimensions, cubes.h(), cubes.mass(), velocity,
                         cubes.mv(), cubes.epsilon(), cubes.sigma(),
                         cubes.type());

      cg.generate(particles);
    }

    for (const auto& spheres : config->spheroids()->spheroid()) {
      const auto& _origin = spheres.origin();
      const auto& _velocity = spheres.velocity();
      dvec3 origin = {_origin.x(), _origin.y(), _origin.z()};
      dvec3 velocity = {_velocity.x(), _velocity.y(), _velocity.z()};
      const int radius = spheres.radius();

      SpheroidGenerator sg(origin, radius, spheres.h(), spheres.mass(),
                           velocity, spheres.epsilon(), spheres.sigma(),
                           spheres.type(), true);

      sg.generate(particles);
    }
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
  }
}
