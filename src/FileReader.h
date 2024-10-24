/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once


#include "defs/Particle.h"

#include <list>

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(std::list<Particle> &particles, const char *filename);
};
