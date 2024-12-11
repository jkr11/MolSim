//
// Created by mcarn on 12/11/24.
//

#pragma once
#include <vector>
#include <memory>

#include "SingularForce.h"
#include "InteractiveForce.h"
#include "ParticleContainer.h"

class VerletIntegrator final {
 protected:
  std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces;
  std::vector<std::unique_ptr<SingularForce>>& singular_forces;
  double delta_t;

 public:
  VerletIntegrator(
      std::vector<std::unique_ptr<InteractiveForce>>& interactive_forces,
      std::vector<std::unique_ptr<SingularForce>>& singular_forces,
      const double delta_t)
      : interactive_forces(interactive_forces),
        singular_forces(singular_forces),
        delta_t(delta_t)
  {}

  ~VerletIntegrator() = default;

  void step(ParticleContainer& particle_container);
};


