//
// Created by jkr on 12/3/24.
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#pragma once
#include "Simulation.h"
#include "containers/ParticleContainer.h"

class Thermostat {
 public:
  double T_init{};
  double T_target{};
  double deltaT{};
  int n_thermostat{};

  explicit Thermostat(const ThermostatConfig& config);

  static double getTemperature(ParticleContainer& particle_container);
  void setTemperature(ParticleContainer& particle_container) const;
};
#endif  // THERMOSTAT_H
