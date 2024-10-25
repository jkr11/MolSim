#pragma once
#include "../defs/Particle.h"
#include "../defs/ParticleContainer.h"
#include "../forces/force.h"
#include "../utils/ArrayUtils.h"
/**
 * @brief Interface for different types of integrators, computing one step in time
 */
class Integrator {
protected:
  Force& force;
  double delta_t;

public:
  Integrator(Force& _force, const double _delta_t) : force{_force}, delta_t{_delta_t} {}
  virtual ~Integrator() {}
  // maybe don't use environment here, I think we only need dt and other things can both be pair or single forces eg Gravity, lennard_jones etc
  virtual void step(ParticleContainer& particle_container) = 0;
};
