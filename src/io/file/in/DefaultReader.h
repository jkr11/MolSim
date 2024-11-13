/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <vector>

#include "defs/Particle.h"
#include "io/file/in/FileReader.h"

/**
 * @brief reads files in the default format, so single line particles, use this
 * for week 1
 */
class DefaultReader final : public FileReader {
 public:
  /**
   * @brief reads files in the default format
   * @param particles vector passed by particle_container
   * @param filename path to the file to be read
   */
  void read(std::vector<Particle>& particles,
                   const std::string& filename) override;
};
