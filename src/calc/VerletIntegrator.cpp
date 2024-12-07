//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
  });

  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  particle_container.singleIterator(
      [this](Particle& p) { p.setF(p.getF() + singularForce.applyForce(p)); });

  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    const dvec3 g12 = interactive_force.directionalForce(p1, p2);
    p1.setF(p1.getF() + g12);  // F_i = \sum_j F_ij
    p2.setF(p2.getF() - g12);  // g12 = -g21
  });

  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
