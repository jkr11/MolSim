//
// Created by jkr on 12/3/24.
//

#include "Thermostat.h"

#include "containers/ParticleContainer.h"
#include "utils/ArrayUtils.h"

Thermostat::Thermostat(const ThermostatConfig &config) {
  this->T_init = config.T_init;
  this->T_target = config.T_target;
  this->deltaT = config.deltaT;
  this->n_thermostat = config.n_thermostat;
}

double Thermostat::getTemperature(ParticleContainer &particle_container) {
  constexpr double D = 2;
  const auto E_kin = particle_container.getKineticEnergy();
  return (2 * E_kin) /
         (D * static_cast<double>(particle_container.getParticles().size()));
}

void Thermostat::setTemperature(ParticleContainer &particle_container) const {
  const double currentTemp = getTemperature(particle_container);
  SpdWrapper::get()->info("Current temp erature is {}", currentTemp);
  const double dT = T_target - currentTemp;
  SpdWrapper::get()->info("Gradient is {} < {}", dT, deltaT);
  if (std::abs(dT) < deltaT) {
    const double new_t = currentTemp + dT;
    SpdWrapper::get()->info("New temperature is {}", new_t);
    double beta = std::sqrt(new_t / currentTemp);
    particle_container.singleIterator(
        [beta](Particle &p) { p.setV(beta * p.getV()); });
  } else {
    const double new_t = currentTemp + deltaT;
    double beta = std::sqrt(new_t / currentTemp);
    particle_container.singleIterator(
        [beta](Particle &p) { p.setV(beta * p.getV()); });
  }
}