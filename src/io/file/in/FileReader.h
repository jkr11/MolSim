#pragma once

#include "defs/Particle.h"

/**
 * @brief abstract class for File Readers of all kinds
 */
class FileReader {
 public:
  /**
   * @brief reads the file into the vector provided by the container
   * @param particles vector passed from particle_container
   * @param filepath path of the file to be read
   */
  static void read(std::vector<Particle> &particles,
                   const std::string &filepath);
};
