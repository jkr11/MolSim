#pragma once
#include <memory>
#include <vector>

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

  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  std::vector<std::unique_ptr<SingularForce>> singular_forces;
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
             const double delta_t,
             std::vector<std::unique_ptr<InteractiveForce>> interactive_forces,
             std::vector<std::unique_ptr<SingularForce>> singular_forces)
      : interactive_force(interactive_force),
        singularForce(singular_force),
        interactive_forces(std::move(interactive_forces)),
        singular_forces(std::move(singular_forces)),
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
