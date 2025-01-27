//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include <iostream>

#include "../utils/ArrayUtils.h"
#include "utils/SpdWrapper.h"

void VerletIntegrator::step(ParticleContainer& particle_container) const {
  // position update
  auto p = particle_container.getParticles()[0];
  std::cout << "Particle " << p->getId() << " is at mem " << p << std::endl;
  auto p2 = particle_container.getParticles()[1];
  for (auto [diag, ref] : p2->getNeighbours()) {
    auto dings = reinterpret_cast<Particle*>(ref);
    if (dings->getId() == 0) {
      std::cout << "Neighbour 0 from Particle 1 has reference location " << ref
                << " and is diag? " << diag << std::endl;
    }
  }

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
  particle_container.computeIndexForces(index_forces_);

  // Lennard Jones (or truncated)
  particle_container.computeInteractiveForces(interactive_forces_);

  // Gravity and or Membrane
  particle_container.computeSingularForces(singular_forces_);

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
