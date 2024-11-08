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

void CuboidReader::read(std::vector<Particle>& particles,
                        const std::string& fileName) {
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
    size_t lineNum = 0;
    std::getline(inputfile, line);
    lineNum++;
    while (line.empty() or line[0] == '#') {
      getline(inputfile, line);
      lineNum++;
      SpdWrapper::get()->debug("Reading line {}: {}", lineNum, line);
    }
    // std::getline(inputfile, line);
    std::istringstream numstream(line);
    numstream >> ncubes;
    SpdWrapper::get()->debug("Reading {} cuboids from file {}", ncubes,
                             fileName);
    if (ncubes == 0) {
      SpdWrapper::get()->error("No cuboids in file {}", fileName);
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < ncubes; i++) {
      std::getline(inputfile, line);
      lineNum++;
      std::istringstream linestream(line);
      // so any of these error dont really make alot of sense because of how we
      // are reading values, if one is missing the reader just skips to the next
      // values and we get an error at the end - missing items. This could be
      // solved by tokenizing our inputfile, but thats too much work for now,
      // mabye this becomes relevant when reading larger files.
      for (auto& c : corner) {
        if (!(linestream >> c)) {
          SpdWrapper::get()->error("Failed to read coordinates in line {}",
                                   lineNum);
          exit(EXIT_FAILURE);
        }
        SpdWrapper::get()->debug("Reading coordinates {}", c);
      }
      // SpdWrapper::get()->debug("corner: {}", corner);

      for (auto& v : velocity) {
        if (!(linestream >> v)) {
          SpdWrapper::get()->error("Failed to read velocity in line {}",
                                   lineNum);
          exit(EXIT_FAILURE);
        }
        SpdWrapper::get()->debug("Reading velocity {}", v);
      }
      // SpdWrapper::get()->debug("velocity {}", velocity);

      for (auto& d : dimensions) {
        if (!(linestream >> d)) {
          SpdWrapper::get()->error("Failed to read dimensions in line {}",
                                   lineNum);
          exit(EXIT_FAILURE);
        }
        SpdWrapper::get()->debug("Reading dimensions {}", d);
      }
      // SpdWrapper::get()->debug("dimensions: {}", dimensions);

      if (!(linestream >> mass)) {
        SpdWrapper::get()->error("Error reading mass at line {}", lineNum);
        exit(EXIT_FAILURE);
      }
      SpdWrapper::get()->debug("mass {} at line {}", mass, lineNum);

      if (!(linestream >> type)) {
        SpdWrapper::get()->error("Error reading type at line {}", lineNum);
        exit(EXIT_FAILURE);
      }

      if (!(linestream >> h)) {
        SpdWrapper::get()->error("Error reading h at line {}", lineNum);
        exit(EXIT_FAILURE);
      }
      SpdWrapper::get()->debug("h: {} at line {}", h, lineNum);

      if (!(linestream >> mv)) {
        SpdWrapper::get()->error("Error reading mv at line {}", lineNum);
        exit(EXIT_FAILURE);
      }
      SpdWrapper::get()->debug("mv: {} at line {}", mv, lineNum);

      if (!(linestream >> epsilon)) {
        SpdWrapper::get()->error("Error reading epsilon at line {}", lineNum);
        exit(EXIT_FAILURE);
      }
      SpdWrapper::get()->debug("epsilon: {} at line {}", epsilon, lineNum);

      if (!(linestream >> sigma)) {
        SpdWrapper::get()->error("Error reading sigma at line {}", lineNum);
        exit(EXIT_FAILURE);
      }
      SpdWrapper::get()->debug("sigma: {} at line {}", sigma, lineNum);

      CuboidGenerator cg(corner, dimensions, h, mass, velocity, mv, epsilon,
                         sigma, type);

      cg.generate(particles);
      SpdWrapper::get()->debug("particle container size {}", particles.size());
    }
  }
}
