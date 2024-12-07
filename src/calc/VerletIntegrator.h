//
// Created by jkr on 10/18/24.
//

#ifndef VERLET_H
#define VERLET_H
#include "Integrator.h"

/**
 * @brief Integrator using the Verlet-integration method on a particle
 * container.
 */
class VerletIntegrator final : public Integrator {
 public:
  /**
   * @brief Create VerletIntegrator object
   * @param interactive_force Reference to the type of force applied each
   * @param singular_force Global force acting on each particle as is
   * @param delta_t Delta time
   */
  VerletIntegrator(InteractiveForce& interactive_force,
                   SingularForce& singular_force, const double delta_t)
      : Integrator(interactive_force, singular_force, delta_t) {}

  /**
   * @brief Destructor
   */
  ~VerletIntegrator() override = default;

  /**
   * @brief Advance particle-system by one time-step.
   * @note Update order: \n
   *   1) Particle positions and boundary conditions \n
   *   2) Intra particular forces \n
   *   3) Particle velocities \n
   *
   * @param particle_container Reference to the current particle-system
   */
  void step(ParticleContainer& particle_container) override;
};
#endif  // VERLET_H
