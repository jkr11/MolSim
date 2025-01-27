//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include <iostream>

#include "../utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

void VerletIntegrator::step(ParticleContainer& particle_container) const {
  // position update
  particle_container.singleIterator([this](Particle& p) {
    // ignore position update for walls
    if (p.getType() < 0) {
      return;
    }

    dvec3 oldx = p.getX();
    const dvec3 new_x = p.getX() + delta_t_ * p.getV() +
                        (delta_t_ * delta_t_ / (2 * p.getM())) * (p.getF());

    p.setX(new_x);
  });

  // reset
  particle_container.singleIterator([](Particle& p) { p.updateForceInTime(); });

  particle_container.imposeInvariant();

  // TODO: refactor in lower iterator? maybe pass time to all? just get the
  // global var?
  // Pull Up
  //particle_container.computeIndexForces(index_forces_);

  particle_container.singleIterator([this](Particle& p) {
  dvec3 f = {0, 0, 0};
  for (const auto& index_force : index_forces_) {
    for (const int id : index_force->getIndeces()) {
      if (p.getId() == id) {
        // SpdWrapper::get()->info("Particle {}; -> [{}, {}, {}]", p.getId(),
        // f[0], f[1], f[2]);
        f = f + index_force->applyForce(p, current_time_);
        // SpdWrapper::get()->info("Particle {} adding [{}, {}, {}]",
        // p.getId(), f[0], f[1], f[2]);
      }
    }
  }
  p.addF(f);
});

  // Lennard Jones (or truncated)cd
  particle_container.computeInteractiveForces(interactive_forces_);

  // Gravity and or Membrane
  //particle_container.computeSingularForces(singular_forces_);
  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& force : singular_forces_) {
      f = f + force->applyForce(p);
    }

    p.addF(f);
  });


  // Velocity Update
  particle_container.singleIterator([this](Particle& p) {
    // ignore position update for walls
    if (p.getType() < 0) {
      return;
    }

    const dvec3 new_v =
        p.getV() + (delta_t_ / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });

  // TODO: remove
  particle_container.incrementTime();
}
