/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <list>

#include "defs/Particle.h"

class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  static void readFile(std::list<Particle>& particles,
                       const std::string& filename);
};
