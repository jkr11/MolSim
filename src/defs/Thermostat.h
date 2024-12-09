//
// Created by jkr on 12/3/24.
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#pragma once
#include "Simulation.h"
#include "containers/ParticleContainer.h"
/**
 * adjusts the temperature of the system at a given periodic frequency
 */
class Thermostat {
 public:
  double T_init{};
  double T_target{};
  double d_temp{};
  int n_thermostat{};
  bool use_relative{};

  explicit Thermostat(const ThermostatConfig& config);

  /**
   * @brief calculates the temperature of the entire system
   * @param particle_container container with the particle system
   * @return the temperature calculated as (2 * E_kin) / 2 * (#dim * #particles)
   */
  static double getTemperature(ParticleContainer& particle_container);

  /**
   * @calculates the average global velocity
   * @param particle_container container
   * @return the average global velocity
   */
  static dvec3 getGlobalVelocity(ParticleContainer& particle_container);

  /**
   * @brief scales the temperature relatively or absolutely using beta
   * @param particle_container container
   * @param beta the scaling coefficient sqrt(T_new / T_current)
   */
  void applyBeta(ParticleContainer& particle_container, double beta) const;

  /**
   * @a wrapper to proved applyBeta with the necessary data
   * @param particle_container container
   */
  void setTemperature(ParticleContainer& particle_container) const;
};
#endif  // THERMOSTAT_H
