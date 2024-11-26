/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include <fstream>
#include <list>

#include "defs/Particle.h"

namespace outputWriter {

class XYZWriter final {
 public:
  XYZWriter();

  ~XYZWriter();

  static void plotParticles(const std::list<Particle>& particles,
                            const std::string& filename, int iteration);
};

}  // namespace outputWriter
