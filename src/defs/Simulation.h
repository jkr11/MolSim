//
// Created by jkr on 11/17/24.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#pragma once
#include <stdexcept>

#include "defs/Particle.h"
#include "utils/SpdWrapper.h"

struct Simulation {
  double delta_t{};
  double t_end{};
  double cutoff_radius{};
  ivec3 domain{};
};

template <typename SVec, typename TVec>
TVec unwrapVec(const SVec& source, const std::string& paramName) {
  try {
    return TVec{source.x(), source.y(), source.z()};
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to unwrap vector " + paramName + ": " +
                             e.what());
  }
}

#endif  // SIMULATION_H
