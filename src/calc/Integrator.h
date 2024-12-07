#pragma once
#include "defs/containers/ParticleContainer.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularForce.h"

/**
 * @brief Interface for different types of integrators
 */
class Integrator {
 protected:
  InteractiveForce &interactive_force;
  SingularForce &singularForce;
  double delta_t;

 public:
  /**
   * @brief Create Integrator object
   * @param interactive_force Reference to the type of force applied each
   * iteration
   * @param singular_force singular force acting on single particles globally
   * @param delta_t Delta time
   * @note Since this is an interface, it's invalid
   */
  Integrator(InteractiveForce &interactive_force, SingularForce &singular_force,
             const double delta_t)
      : interactive_force(interactive_force),
        singularForce(singular_force),
        delta_t(delta_t) {};

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
