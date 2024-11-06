//
// Created by jkr on 11/1/24.
//
#include "CuboidReader.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "defs/CuboidGenerator.h"
#include "defs/Particle.h"
#include "utils/SpdWrapper.h"

void CuboidReader::readCuboidFile(ParticleContainer& particle_container,
                                  std::string& fileName) {
  std::ifstream inputfile(fileName);
  SpdWrapper::get()->debug("Reading cuboid file {}", fileName);
  if (inputfile.is_open()) {
    dvec3 corner{};
    dvec3 velocity{};
    std::array<int, 3> dimensions{};
    int ncubes;
    int type;
    double h;
    double mass;
    double mv;
    double epsilon;
    double sigma;
    std::string line;
    std::getline(inputfile, line);
    while (line.empty() or line[0] == '#') {
      getline(inputfile, line);
      SpdWrapper::get()->debug("Reading line {}", line);
    }
    // std::getline(inputfile, line);
    std::istringstream numstream(line);
    numstream >> ncubes;
    SpdWrapper::get()->debug("Reading {} cuboids from file {}", ncubes,
                             fileName);
    if (ncubes == 0) {
      SpdWrapper::get()->error("No cuboids in file {}", fileName);
    }
    for (int i = 0; i < ncubes; i++) {
      std::getline(inputfile, line);
      std::istringstream linestream(line);
      for (auto& c : corner) {
        linestream >> c;
      }
      for (auto& v : velocity) {
        linestream >> v;
      }
      for (auto& d : dimensions) {
        linestream >> d;
      }
      if (linestream.eof()) {
        SpdWrapper::get()->error("Error reading file: eof reached");
        exit(-1);
      }
      linestream >> mass;
      SpdWrapper::get()->debug("mass {}", mass);
      linestream >> type;
      SpdWrapper::get()->debug("type: {}", type);
      linestream >> h;
      SpdWrapper::get()->debug("h: {}", h);
      linestream >> mv;
      SpdWrapper::get()->debug("mv: {}", mv);
      linestream >> epsilon;
      SpdWrapper::get()->debug("epsilon: {}", epsilon);
      linestream >> sigma;
      SpdWrapper::get()->debug("sigma: {}", sigma);
      CuboidGenerator cg(corner, dimensions, h, mass, velocity, mv, epsilon,
                         sigma, type);
      cg.generate(particle_container);
      SpdWrapper::get()->debug("particle container size {}",
                               particle_container.size());
    }
  }
}
