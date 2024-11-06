/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <vector>
#include "defs/Particle.h"
#include "io/file/FileReader.h"

class DefaultReader : FileReader {
 public:
  static void read(std::vector<Particle>& particles,
                       const std::string& filename);
};
