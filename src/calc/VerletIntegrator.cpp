//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include <iostream>

#include "../utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      return;
    }  // ignore position update for walls
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
  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& index_force : index_forces_) {
      for (const int id : index_force->getIndeces()) {
        if (p.getId() == id) {
          f = f + index_force->applyForce(p, current_time_);
        }
      }
    }
    p.addF(f);
  });

  // Lennard Jones (or truncated)
  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    dvec3 f12 = {0.0, 0.0, 0.0};
    for (const auto& force : interactive_forces_) {
      f12 = f12 + force->directionalForce(p1, p2);
    }
    p1.addF(f12);
    p2.subF(f12);
  });

  // Gravity and or Membrane
  particle_container.singleIterator([this](Particle& p) {
    dvec3 f = {0, 0, 0};
    for (const auto& force : singular_forces_) {
      f = f + force->applyForce(p);
    }
    p.addF(f);
  });

  // Velocity Update
  particle_container.singleIterator([this](Particle& p) {
    if (p.getType() < 0) {
      return;  // ignore velocity update for walls, theoretically
    }

    const dvec3 new_v =
        p.getV() + (delta_t_ / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });

  // TODO: move to container (specify const again)
  incrementTime();
}
