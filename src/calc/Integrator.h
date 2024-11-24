#pragma once
#include "../defs/containers/ParticleContainer.h"
#include "../forces/BidirectionalForce.h"

/**
 * @brief Interface for different types of integrators
 */
class Integrator {
 protected:
  BidirectionalForce &force;
  double delta_t;

 public:
  /**
   * @brief Create Integrator object
   * @param force Reference to the type of force applied each iteration
   * @param delta_t Delta time
   * @note Since this is an interface, it's invalid
   */
  Integrator(BidirectionalForce &force, const double delta_t)
      : force(force), delta_t(delta_t) {};

  /**
   * @brief Virtual destructor for all Integrator inheritors
   */
  virtual ~Integrator() = default;

  /**
   * @brief Virtual method to advance time by one step
   * @param particle_container
   */
  virtual void step(ParticleContainer &particle_container) = 0;
};
