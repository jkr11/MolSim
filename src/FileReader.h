/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

<<<<<<< HEAD
#include "Particle.h"
=======
#include "defs/Particle.h"
>>>>>>> MolSim/test1

#include <list>

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  void readFile(std::list<Particle> &particles, char *filename);
};
