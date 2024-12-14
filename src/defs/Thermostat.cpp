//
// Created by jkr on 12/3/24.
//

#include "Thermostat.h"

#include "containers/ParticleContainer.h"
#include "utils/ArrayUtils.h"

Thermostat::Thermostat(const ThermostatConfig &config) {
  this->T_init = config.T_init;
  this->T_target = config.T_target;
  this->d_temp = config.deltaT;
  this->n_thermostat = config.n_thermostat;
  this->use_relative = config.use_relative;
}

double Thermostat::getTemperature(ParticleContainer &particle_container) {
  constexpr double D = 2;
  const auto E_kin = particle_container.getKineticEnergy();
  return (2 * E_kin) /
         (D * static_cast<double>(particle_container.getParticles().size()));
}

dvec3 Thermostat::getGlobalVelocity(ParticleContainer &particle_container) {
  dvec3 avg_velocity_acc = {0.0, 0.0, 0.0};
  particle_container.singleIterator([&avg_velocity_acc](const Particle &p) {
    avg_velocity_acc = avg_velocity_acc + p.getV();
  });
  return avg_velocity_acc;
}
void Thermostat::applyBeta(ParticleContainer &particle_container,
                           const double beta) const {
  const auto average_velocity = getGlobalVelocity(particle_container);
  if (use_relative) {
    particle_container.singleIterator([&average_velocity, &beta](Particle &p) {
      p.subV(average_velocity);
      p.mulV(beta);
      p.addV(average_velocity);
    });
  } else {
    particle_container.singleIterator([&beta](Particle &p) { p.mulV(beta); });
  }
}

void Thermostat::setTemperature(ParticleContainer &particle_container) const {
  const double current_temp = getTemperature(particle_container);
  SpdWrapper::get()->info("current temperature is {}", current_temp);
  const double dT = T_target - current_temp;

  const auto adjustment = std::abs(dT) < d_temp ? dT : d_temp;
  const double new_temp = current_temp + adjustment;
  SpdWrapper::get()->info("new_temp is {}", new_temp);
  const double beta = std::sqrt(new_temp / current_temp);
  SpdWrapper::get()->info("beta is {}", beta);
  applyBeta(particle_container, beta);
}