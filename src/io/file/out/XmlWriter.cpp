//
// Created by jkr on 12/3/24.
//
#include "XmlWriter.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "defs/Thermostat.h"
#include "io/file/out/checkpoint-schema.hxx"

XmlWriter::XmlWriter() = default;

XmlWriter::~XmlWriter() = default;

template <typename VecType, typename XmlType>
inline XmlType wrapVec(const VecType& vec) {
  XmlType xmlVec{vec[0], vec[1], vec[2]};
  return xmlVec;
}

ParticleType wrapParticle(const Particle& particle) {
  const auto position = wrapVec<dvec3, CDvec3Type>(particle.getX());
  const auto velocity = wrapVec<dvec3, CDvec3Type>(particle.getV());
  const auto force =
      wrapVec<std::array<double, 3>, CDvec3Type>(particle.getF());
  const auto old_force =
      wrapVec<std::array<double, 3>, CDvec3Type>(particle.getOldF());
  const auto mass = particle.getM();
  const auto epsilon = particle.getEpsilon();
  const auto sigma = particle.getSigma();
  const auto type = particle.getType();
  return ParticleType{position, velocity, force, old_force,
                      mass,     epsilon,  sigma, type};
}

void XmlWriter::writeFile(ParticleContainer& particle_container,
                          const std::string& filepath) {
  try {
    ParticlesType xml_particles;

    particle_container.singleIterator(
        [&xml_particles](const Particle& particle) {
          const ParticleType xml_particle = wrapParticle(particle);
          xml_particles.Particle().push_back(xml_particle);
        });
    xml_schema::namespace_infomap map;

    map[""].name = "";
    map[""].schema = "../../src/io/file/out/checkpoint-schema.xsd";
    // I think this can either stay hardcoded or we provide a symlink

    std::ostringstream fileName;
    CheckpointType checkpoint{xml_particles};

    std::ofstream checkpoint_file(filepath);
    SpdWrapper::get()->info("--- Written checkpoint to {}", filepath);
    Checkpoint(checkpoint_file, checkpoint, map);
    checkpoint_file.close();

  } catch (const std::exception& e) {
    std::cerr << "Error writing XML file: " << e.what() << "\n";
  }
}