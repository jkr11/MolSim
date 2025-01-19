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
  this->use_thermal_motion = config.use_thermal_motion;
}

double Thermostat::getTemperature(ParticleContainer &particle_container) {
  constexpr double D = 2;  // TODO: make global so its 3D
  const auto E_kin = particle_container.getKineticEnergy();
  return (2 * E_kin) /
         (D * static_cast<double>(particle_container.getParticles().size()));
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
                           const double beta) const {
  particle_container.singleIterator([&beta](Particle &p) { p.mulV(beta); });
}

void Thermostat::applyThermalBeta(ParticleContainer &particle_container,
                                  const double beta,
                                  const dvec3 &avg_velocity) const {
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

  return {total_velocity[0] / c, total_velocity[1] / c, total_velocity[2 / c]};
}

double Thermostat::getThermalTemperature(ParticleContainer &particle_container,
                                         dvec3 avg_velocity) {
  constexpr double D = 3;
  double E_kin = 0.0;
  particle_container.singleIterator([&E_kin, avg_velocity](const Particle &p) {
    E_kin += p.getM() * ArrayUtils::L2InnerProduct(p.getV() - avg_velocity);
  });

  return E_kin / (D * static_cast<double>(
                          particle_container.getParticleCount() -
                          particle_container.getSpecialParticleCount()));
}

void Thermostat::setTemperature(ParticleContainer &particle_container) const {
  if (use_thermal_motion) {
    const dvec3 average_velocity = getAverageVelocity(particle_container);
    const double current_temp =
        getThermalTemperature(particle_container, average_velocity);
#ifndef BENCHMARK
    SpdWrapper::get()->info("current temperature is {}", current_temp);
#endif
    const double dT = T_target - current_temp;

    double adjustment = 0;
    if (std::abs(dT) > d_temp) {
      adjustment = (dT / std::abs(dT)) * d_temp;
    } else {
      adjustment = dT;
    }
#ifndef BENCHMARK
    SpdWrapper::get()->info("adjustment is {}", adjustment);
#endif
    const double new_temp = current_temp + adjustment;
#ifndef BENCHMARK
    SpdWrapper::get()->info("new_temp is {}", new_temp);
#endif
    const double beta = std::sqrt(new_temp / current_temp);
#ifndef BENCHMARK
    SpdWrapper::get()->info("beta is {}", beta);
#endif
    applyThermalBeta(particle_container, beta, average_velocity);
#ifndef BENCHMARK
    auto temp_after = getTemperature(particle_container);
    SpdWrapper::get()->info("temp_after is {}", temp_after);
#endif

  } else {
    const double current_temp = getTemperature(particle_container);
#ifndef BENCHMARK
    SpdWrapper::get()->info("current temperature is {}", current_temp);
#endif
    const double dT = T_target - current_temp;

    double adjustment = 0;
    if (std::abs(dT) > d_temp) {
      adjustment = (dT / std::abs(dT)) * d_temp;
    } else {
      adjustment = dT;
    }
#ifndef BENCHMARK
    SpdWrapper::get()->info("adjustment is {}", adjustment);
#endif
    const double new_temp = current_temp + adjustment;
#ifndef BENCHMARK
    SpdWrapper::get()->info("new_temp is {}", new_temp);
#endif
    const double beta = std::sqrt(new_temp / current_temp);
#ifndef BENCHMARK
    SpdWrapper::get()->info("beta is {}", beta);
#endif
    applyBeta(particle_container, beta);
#ifndef BENCHMARK
    auto temp_after = getTemperature(particle_container);
    SpdWrapper::get()->info("temp_after is {}", temp_after);
#endif
  }
}