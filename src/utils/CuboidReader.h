//
// Created by jkr on 11/1/24.
//
#pragma once
#ifndef CUBOIDREADER_H
#define CUBOIDREADER_H
#include <string>
#include "defs/ParticleContainer.h"
class CuboidReader {
  public:
  CuboidReader() = default;
  ~CuboidReader() = default;
  static void readCuboidFile(ParticleContainer& particle_container, std::string& fileName);
};
#endif //CUBOIDREADER_H
