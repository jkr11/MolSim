//
// Created by jkr on 10/18/24.
//
#include "VerletIntegrator.h"

#include <iostream>

#include "../utils/ArrayUtils.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  // update positions
  particle_container.single_iterator([this](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() +
                        (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
  });

  // advance time, old F_{t] -> F_{t-1}, F_{t} = 0
  particle_container.single_iterator([](Particle& p) {
    p.setOldF(p.getF());
    p.setF({0, 0, 0});
  });

  // update force
  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    const dvec3 g12 = force.directionalForce(p1, p2);
    p1.setF(p1.getF() + g12);  // F_i = \sum_j F_ij
    p2.setF(p2.getF() - g12);  // g12 = -g21
  });
  // Now we use F_t and F_{t-1} to calculate the current velocity
  particle_container.single_iterator([this](Particle& p) {
    const dvec3 new_v =
        p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
