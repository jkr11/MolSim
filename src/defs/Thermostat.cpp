//
// Created by jkr on 12/3/24.
//

#include "Thermostat.h"

#include "containers/ParticleContainer.h"
#include "debug/debug_print.h"
#include "utils/ArrayUtils.h"

Thermostat::Thermostat(const ThermostatConfig &config) {
  this->t_init_ = config.t_init;
  this->t_target_ = config.t_target;
  this->delta_temp_ = config.delta_t;
  this->n_thermostat_ = config.n_thermostat;
  this->use_relative_ = config.use_relative;
  this->use_thermal_motion_ = config.use_thermal_motion;
  this->dimension_ = config.two_d ? 2 : 3;
}

double Thermostat::getTemperature(ParticleContainer &particle_container) const {
  const auto e_kin = particle_container.getKineticEnergy();
  return (2 * e_kin) /
         (dimension_ *
          static_cast<double>(particle_container.getParticles().size()));
}

dvec3 Thermostat::getGlobalVelocity(ParticleContainer &particle_container) {
  dvec3 avg_velocity_acc = {0.0, 0.0, 0.0};
  particle_container.singleIterator([&avg_velocity_acc](const Particle &p) {
    avg_velocity_acc = avg_velocity_acc + p.getV();
  });
  return (1.0 / static_cast<double>(particle_container.size())) *
         avg_velocity_acc;
}

void Thermostat::applyBeta(ParticleContainer &particle_container,
                           const double beta) {
  particle_container.singleIterator([&beta](Particle &p) { p.mulV(beta); });
}

void Thermostat::applyThermalBeta(ParticleContainer &particle_container,
                                  const double beta,
                                  const dvec3 &avg_velocity) {
  particle_container.singleIterator([&beta, avg_velocity](Particle &p) {
    p.setV(avg_velocity + beta * (p.getV() - avg_velocity));
  });
}

dvec3 Thermostat::getAverageVelocity(ParticleContainer &particle_container) {
  dvec3 total_velocity = {0.0, 0.0, 0.0};
  const auto c =
      static_cast<double>(particle_container.getParticleCount() -
                          particle_container.getSpecialParticleCount());

  for (const auto &p : particle_container.getParticles()) {
    if (p->getType() >= 0) {
      total_velocity = total_velocity + p->getV();
    }
  }

  return {total_velocity[0] / c, total_velocity[1] / c, total_velocity[2] / c};
}

double Thermostat::getThermalTemperature(ParticleContainer &particle_container,
                                         dvec3 avg_velocity) const {
  double e_kin = 0.0;
  particle_container.singleIterator([&e_kin, avg_velocity](const Particle &p) {
    if (p.getType() < 0) return;  // exclude walls
    e_kin += p.getM() * ArrayUtils::L2InnerProduct(p.getV() - avg_velocity);
  });

  return e_kin /
         (dimension_ *
          static_cast<double>(particle_container.getParticleCount() -
                              particle_container.getSpecialParticleCount()));
}

void Thermostat::setTemperature(ParticleContainer &particle_container) const {
  dvec3 average_velocity;
  double current_temp = 0;

  if (use_thermal_motion_) {
    average_velocity = getAverageVelocity(particle_container);
    current_temp = getThermalTemperature(particle_container, average_velocity);
  } else {
    current_temp = getTemperature(particle_container);
  }
  DEBUG_PRINT_FMT("current temperature is {}", current_temp);
  const double d_t = t_target_ - current_temp;

  double adjustment = 0;
  if (std::abs(d_t) > delta_temp_) {
    adjustment = (d_t / std::abs(d_t)) * delta_temp_;
  } else {
    adjustment = d_t;
  }
  DEBUG_PRINT_FMT("adjustment is {}", adjustment);
  const double new_temp = current_temp + adjustment;
  DEBUG_PRINT_FMT("new_temp is {}", new_temp);
  const double beta = std::sqrt(new_temp / current_temp);
  DEBUG_PRINT_FMT("beta is {}", beta);
  if (use_thermal_motion_) {
    applyThermalBeta(particle_container, beta, average_velocity);
  } else {
    applyBeta(particle_container, beta);
  }
#ifndef BENCHMARK
#ifdef DEBUG
  average_velocity = getAverageVelocity(particle_container);
  auto temp_after = getThermalTemperature(particle_container, average_velocity);
  DEBUG_PRINT_FMT("temp_after is {}", temp_after);
#endif
#endif
}

int Thermostat::getNThermostat() const { return n_thermostat_; }

double Thermostat::getTTarget() const { return t_target_; }
