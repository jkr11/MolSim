/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

<<<<<<< HEAD
#include "Particle.h"
=======
#include "defs/Particle.h"
>>>>>>> MolSim/test1

#include <fstream>
#include <list>

namespace outputWriter {

class XYZWriter {

public:
  XYZWriter();

  virtual ~XYZWriter();

  void plotParticles(std::list<Particle> particles, const std::string &filename,
                     int iteration);
};

} // namespace outputWriter
