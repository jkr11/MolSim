#include "../defs/Particle.h"
#include "../defs/ParticleContainer.h"
#include "../forces/Force.h"
#include "../utils/ArrayUtils.h"
#pragma once

/**
 * @brief Interface for different types of integrators
 */
class Integrator {
 protected:
  Force& force;
  double delta_t;

 public:
  /**
   * @brief Create Integrator object
   * @param force Reference to the type of force applied each iteration
   * @param delta_t Delta time
   * @note Since this is an interface, it's invalid
   */
  Integrator(Force& force, double delta_t) : force(force), delta_t(delta_t) {};

  /**
   * @brief Virtual destructor for all Integrator inheritors
   */
  virtual ~Integrator() = default;

  /**
   * @brief Virtual method to advance time by one step
   * @param particle_container
   */
  virtual void step(ParticleContainer& particle_container) = 0;
};
