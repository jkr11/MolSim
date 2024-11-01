//
// Created by jkr on 11/1/24.
//
#include <iostream>
#include<sstream>
#include<fstream>
#include "CuboidReader.h"
#include "../defs/Particle.h"
#include "../defs/CuboidGenerator.h"
void CuboidReader::readCuboidFile(ParticleContainer & particle_container, std::string& fileName) {
  std::ifstream inputfile(fileName);
  dvec3 corner{};
  dvec3 velocity{};
  double temp;
  double mass;
  double h;
  int type;
  std::array<int,3> dimensions{};
  if(inputfile.is_open()) {
    std::string line;
    std::getline(inputfile, line);
    std::istringstream linestream(line);
    for(auto& c : corner) {
      linestream >> c;
    }
    for (auto&v : velocity) {
      linestream >> v;
    }
    for (auto& d : dimensions) {
      linestream >> d;
    }
    linestream >> mass;
    std::cout << "mass: " << mass << std::endl;
    linestream >> type;
    std::cout << "type: " << type << std::endl;
    linestream >> h;
    std::cout << "h: " << h << std::endl;
    linestream >> temp;
    std::cout << "temp: " << temp << std::endl;
    CuboidGenerator cg(corner, dimensions,h, mass, velocity, temp, type);
    cg.generate(particle_container);
    std::cout << particle_container.size() << std::endl;
  }
}
