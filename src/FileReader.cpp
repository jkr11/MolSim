/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "debug/debug_print.h"
#include "utils/SpdWrapper.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(std::list<Particle> &particles,
                          const std::string &filename) {
  if (std::ifstream input_file(filename); input_file.is_open()) {
    std::string tmp_string;
    int num_particles = 0;
    double m;
    std::array<double, 3> v{};
    std::array<double, 3> x{};
    int type;

    getline(input_file, tmp_string);
    SpdWrapper::get()->debug("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      SpdWrapper::get()->debug("Read line: {}", tmp_string);
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    SpdWrapper::get()->debug("Reading {}", num_particles);
    getline(input_file, tmp_string);
    SpdWrapper::get()->debug("Read line: {}", tmp_string);

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        SpdWrapper::get()->error("Error reading file: eof reached");
        exit(-1);
      }
      datastream >> m;
      datastream >> type;
      particles.emplace_back(x, v, m, 1.0, 1.0, type);

      getline(input_file, tmp_string);
      SpdWrapper::get()->debug("Read line: {}", tmp_string);
    }
  } else {
    SpdWrapper::get()->error("Error opening file {}", filename);
    exit(-1);
  }
}
