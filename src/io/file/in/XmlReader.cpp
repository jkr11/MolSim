//
// Created by jkr on 11/16/24.
//

#include "XmlReader.h"

#include <fstream>

#include "defs/CuboidGenerator.h"
#include "input.hxx"
#include "utils/SpdWrapper.h"

void XmlReader::read(std::vector<Particle>& particles,
                     const std::string& filepath) {
  try {
    SpdWrapper::get()->info("XMLREADER entered");
    std::ifstream xmlFile(filepath);
    SpdWrapper::get()->info("XMLREADER opened");
    if (!xmlFile.is_open()) {
      SpdWrapper::get()->error("Unable to open file {}", filepath);
    }

    std::unique_ptr<::simulation> config = simulation_(xmlFile);
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
                         static_cast<int>(cubes.type()));

      cg.generate(particles);
    }
  } catch (const std::exception& e) {
    SpdWrapper::get()->error("Error reading XML file: {}", e.what());
  }
}
