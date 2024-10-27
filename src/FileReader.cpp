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

    getline(input_file, tmp_string);
    DEBUG_PRINT("Read line: " + tmp_string + "\n");

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      DEBUG_PRINT("Read line: " + tmp_string + "\n");
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    DEBUG_PRINT("Reading " + num_particles + "." + "\n");
    getline(input_file, tmp_string);
    DEBUG_PRINT("Read line: " + tmp_string + "\n");

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        std::cout
            << "Error reading file: eof reached unexpectedly reading from line "
            << i << std::endl;
        exit(-1);
      }
      datastream >> m;
      particles.emplace_back(x, v, m);

      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}
