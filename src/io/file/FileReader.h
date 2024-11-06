#pragma once

#include <cstdlib>
#include <list>

#include "defs/Particle.h"

class FileReader {
 public:
  static void read(std::vector<Particle>& particles, const std::string& filepath);
};