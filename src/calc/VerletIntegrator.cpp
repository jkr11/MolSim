//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include "../utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& particle_container) const {
  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
  });

  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& force : singular_forces) {
      f = f + force->applyForce(p);
    }
    p.addF(f);
  });

  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    dvec3 f12 = {0.0, 0.0, 0.0};
    for (const auto& force : interactive_forces) {
      f12 = f12 + force->directionalForce(p1, p2);
    }

    //aquire in same order => prevent deadlock (TODO: ?)
    p1.addF(f12);
    p2.subF(f12);
  });

  particle_container.singleIterator([this](Particle& p) {
    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
