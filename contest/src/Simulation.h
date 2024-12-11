//
// Created by mcarn on 12/11/24.
//

#pragma once
#include <vector>
#include <string>

#include "../../src/defs/types.h"

struct LinkedCellsConfig {
  ivec3 domain;
  double cutoff_radius;
  enum BoundaryType { Outflow, Reflective, Periodic } boundary_type;
  struct BoundaryConfig {
    BoundaryType x_high;
    BoundaryType x_low;
    BoundaryType y_high;
    BoundaryType y_low;
    BoundaryType z_high;
    BoundaryType z_low;
  } boundary_config;
};



struct ParticleTypeInfo {
  std::vector<double> mass;
  std::vector<double> sigma;
  std::vector<double> epsilon;
};

class Simulation {
  public:
   static void run(std::string filepath);
};
