//
// Created by jkr on 10/18/24.
//

#ifndef VERLET_H
#define VERLET_H
#pragma once
#include <utility>

#include "defs/containers/ParticleContainer.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularForce.h"
/**
 * @brief Integrator using the Verlet-integration method on a particle
 * container.
 */
class VerletIntegrator {
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces;
  std::vector<std::unique_ptr<SingularForce>> singular_forces;
  double delta_t;

 public:
  /**
   * @brief Create VerletIntegrator object
   * @param interactive_forces Reference to the type of force applied each
   * @param singular_forces Global forces acting on each particle as is
   * @param delta_t Delta time
   */
  VerletIntegrator(
      std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
      std::vector<std::unique_ptr<SingularForce>>& singular_forces,
      const double delta_t)
      : interactive_forces(std::move(interactive_forces)),
        singular_forces(std::move(singular_forces)),
        delta_t(delta_t) {};

  /**
   * @brief Destructor
   */
  ~VerletIntegrator() = default;

  /**
   * @brief Advance particle-system by one time-step.
   * @note Update order: \n
   *   1) Particle positions and boundary conditions \n
   *   2) Intra particular forces \n
   *   3) Particle velocities \n
   *
   * @param particle_container Reference to the current particle-system
   */
  void step(ParticleContainer& particle_container) const;
};
#endif  // VERLET_H
