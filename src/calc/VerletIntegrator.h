#pragma once
#include <memory>
#include <vector>

#include "defs/containers/ParticleContainer.h"
#include "forces/IndexForce.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularForce.h"

/**
 * @brief Interface for different types of integrators
 */
class VerletIntegrator {
 protected:
  std::vector<std::unique_ptr<InteractiveForce>> interactive_forces_;
  std::vector<std::unique_ptr<SingularForce>> singular_forces_;
  std::vector<std::unique_ptr<IndexForce>> index_forces_;
  double delta_t_;
  double current_time_;


 public:
  /**
   * @brief Create Integrator object
   * @param interactive_forces References to the type of force applied each
   * iteration
   * @param singular_forces singular forces acting on single particles globally
   * @param delta_t Delta time
   * @param index_forces
   * @note Since this is an interface, it's invalid
   */
  VerletIntegrator(std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
             std::vector<std::unique_ptr<SingularForce>>& singular_forces,
             const double delta_t,
             std::vector<std::unique_ptr<IndexForce>>& index_forces)
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
   void step(ParticleContainer& particle_container) const;
};
