//
// Created by jkr on 10/18/24.
//

#ifndef VERLET_H
#define VERLET_H
#include "Integrator.h"

class VerletIntegrator : public Integrator {
public:
  VerletIntegrator(Force &force, const double delta_t) : Integrator(force, delta_t) {}
  ~VerletIntegrator() override = default;
  void step(ParticleContainer& particle_container) override;
};
#endif //VERLET_H
