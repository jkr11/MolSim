//
// Created by jkr on 10/18/24.
//

#ifndef VERLET_H
#define VERLET_H
#pragma once

#include "Integrator.h"

/**
 * @brief Integrator using the Verlet-integration method on a particle
 * container.
 */
class VerletIntegrator final : public Integrator {
 public:
  /**
   * @brief Create VerletIntegrator object
   * @param interactive_forces Reference to the type of force applied each
   * @param singular_forces Global forces acting on each particle as is
   * @param index_forces
   * @param delta_t Delta time
   */
  VerletIntegrator(
      std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
      std::vector<std::unique_ptr<SingularForce>>& singular_forces,
      std::vector<std::unique_ptr<IndexForce>>& index_forces,
      const double delta_t)
      : Integrator(interactive_forces, singular_forces, delta_t, index_forces) {
  }

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
