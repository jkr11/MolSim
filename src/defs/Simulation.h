//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once
#include <array>

struct Simulation {
  double delta_t{};
  double t_end{};
  double cutoff_radius{};
  std::array<int, 3> domain{};
};

#endif  // SIMULATION_H
