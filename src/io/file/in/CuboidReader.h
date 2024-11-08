//
// Created by jkr on 11/1/24.
//
#pragma once
#ifndef CUBOIDREADER_H
#define CUBOIDREADER_H

#include <string>

#include "defs/ParticleContainer.h"
#include "io/file/in/FileReader.h"

/**
 * @brief Reads files describing cuboids (.txt, .cuboid) in the following format
 * [TODO: format]
 */
class CuboidReader : FileReader {
 public:
  /**
   * @brief reads the file description of a cuboid into its vector form
   * @param particles is the reference to the vector of particles from a \n
   * particle container, storing the cuboid
   * @param fileName location of the cuboid file
   */
  static void read(std::vector<Particle>& particles,
                   const std::string& fileName);
};
#endif  // CUBOIDREADER_H
