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

  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    // TODO: this should be correct, as there is no cutoff on singular forces
    for (const auto& force : singular_forces) {
      f = f + force->applyForce(p);
    }
    p.setF(p.getF() + f);
  });

  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    dvec3 f12 = {0.0, 0.0, 0.0};
    // TODO: this might either be wrong or inefficient for week5
    for (const auto& force : interactive_forces) {
      f12 = f12 + force->directionalForce(p1, p2);
    }
    p1.setF(p1.getF() + f12);  // F_i = \sum_j F_ij
    p2.setF(p2.getF() - f12);  // g12 = -g21
  });

  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
