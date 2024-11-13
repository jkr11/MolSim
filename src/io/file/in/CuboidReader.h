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
 * [number of cuboids : int]
 * [corner, velocity, dimensions, mass, type, h, mean_velocity, epsilon, sigma]
 */
class CuboidReader final : public FileReader {
 public:
  /**
   * @brief reads the file description of a cuboid into its vector form
   * @param particles is the reference to the vector of particles from a \n
   * particle container, storing the cuboid
   * @param fileName location of the cuboid file
   */
  void read(std::vector<Particle>& particles,
                   const std::string& fileName) override;
};
#endif  // CUBOIDREADER_H
