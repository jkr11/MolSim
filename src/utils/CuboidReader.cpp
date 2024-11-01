//
// Created by jkr on 11/1/24.
//
#include "CuboidReader.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "../defs/CuboidGenerator.h"
#include "../defs/Particle.h"
void CuboidReader::readCuboidFile(ParticleContainer& particle_container,
                                  std::string& fileName) {
  std::ifstream inputfile(fileName);
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
      std::cout << "Read line: " + line << std::endl;
    }
    //sstd::getline(inputfile, line);
    std::istringstream numstream(line);
    numstream >> ncubes;
    std::cout << "Reading " << ncubes << " cuboids from file " << fileName
              << std::endl;
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
      linestream >> mass;
      std::cout << "mass: " << mass << std::endl;
      linestream >> type;
      std::cout << "type: " << type << std::endl;
      linestream >> h;
      std::cout << "h: " << h << std::endl;
      linestream >> mv;
      std::cout << "mv: " << mv << std::endl;
      linestream >> epsilon;
      std::cout << "epsilon: " << epsilon << std::endl;
      linestream >> sigma;
      std::cout << "sigma: " << sigma << std::endl;
      CuboidGenerator cg(corner, dimensions, h, mass, velocity, mv, epsilon,
                         sigma, type);
      cg.generate(particle_container);
      std::cout << particle_container.size() << std::endl;
    }
  }
}
