//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  particle_container.singleIterator(
      [this](Particle& p) { p.updateX(delta_t); });

  /*
  for (Particle& p : particle_container.getParticlesObjects()) {
    const double delta_t_sq_over_2m = delta_t * delta_t / (2 * p.getM());
    const dvec3 new_x =
        p.getX() + delta_t * p.getV() + delta_t_sq_over_2m * (p.getF());
    p.setX(new_x);
    p.updateX(delta_t);
  }
  */

  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    // TODO: this should be correct, as there is no cutoff on singular forces
    for (const auto& force : singular_forces) {
      f = f + force->applyForce(p);
    }
    p.addF(f);
  });

  particle_container.pairIterator(interactive_forces);

  particle_container.singleIterator(
      [this](Particle& p) { p.updateV(delta_t); });
  /*
  for (Particle& p : particle_container.getParticlesObjects()) {
    const double delta_t_sq_over_2m = delta_t * delta_t / (2 * p.getM());
    const dvec3 new_x =
        p.getX() + delta_t * p.getV() + delta_t_sq_over_2m * (p.getF());
    p.setX(new_x);
    p.updateV(delta_t);
  }
*/
}
