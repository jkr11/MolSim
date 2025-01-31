#pragma once

#ifndef VERLETINTEGRATOR_H
#define VERLETINTEGRATOR_H

#include <memory>
#include <vector>

#include "defs/Simulation.h"
#include "forces/IndexForce.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularForce.h"

/**
 * @brief Interface for different types of integrators
 */
class VerletIntegrator {
 protected:
  /**
   * Vector of applied interactive forces
   */
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces_;
  /**
   * Vector of applied singular forces
   */
  std::vector<std::unique_ptr<SingularForce>> singular_forces_;
  /**
   * Vector of applied index forces
   */
  std::vector<std::unique_ptr<IndexForce>> index_forces_;
  /**
   * Time step
   */
  double delta_t_;
  /**
   * Current time for index force
   */
  double current_time_;
  /**
   * Parallelisation strategy used, default is 3 so no parallelisation
   */
  ParallelStrategy strategy_ = ParallelStrategy::STRATEGY_3;

 public:
  /**
   * @brief Create Integrator object
   * @param interactive_forces References to the type of force applied each
   * iteration
   * @param singular_forces singular forces acting on single particles globally
   * @param delta_t Delta time
   * @param strategy whether to use C18 or force buffer
   * @param index_forces Index forces applied each step
   * @note Since this is an interface, it's invalid
   */
  VerletIntegrator(
      std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
      std::vector<std::unique_ptr<SingularForce>>& singular_forces,
      std::vector<std::unique_ptr<IndexForce>>& index_forces,
      const double delta_t, const ParallelStrategy strategy)
      : interactive_forces_(std::move(interactive_forces)),
        singular_forces_(std::move(singular_forces)),
        index_forces_(std::move(index_forces)),
        delta_t_(delta_t),
        current_time_(0),
        strategy_(strategy) {};

  /**
   * @brief Create Integrator object with Strategy 3 (no OpenMP)
   * @param interactive_forces References to the type of force applied each
   * iteration
   * @param singular_forces singular forces acting on single particles globally
   * @param index_forces Index forces applied each step
   * @param delta_t Delta time
   */
  explicit VerletIntegrator(
      std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
      std::vector<std::unique_ptr<SingularForce>>& singular_forces,
      std::vector<std::unique_ptr<IndexForce>>& index_forces,
      const double delta_t)
      : interactive_forces_(std::move(interactive_forces)),
        singular_forces_(std::move(singular_forces)),
        index_forces_(std::move(index_forces)),
        delta_t_(delta_t),
        current_time_(0) {};

  /**
   * @brief Virtual destructor for all Integrator inheritors
   */
  ~VerletIntegrator() = default;

  /**
   * @brief Virtual method to advance time by one step
   * @param particle_container
   */
  void step(ParticleContainer& particle_container);

  /**
   * @brief increments the current simulation time by delta_t_
   */
  void incrementTime() { current_time_ += delta_t_; }
};
#endif