//
// Created by jkr on 10/18/24.
//
#include "Verlet.h"

void VerletIntegrator::step(ParticleContainer& particle_container) {
  //update X -> F -> V
  //[this] needs to be captured everywhere since delta_t is a class constant
  particle_container.single_iterator([this](Particle& p) {
    const dvec3 new_x = p.getX() + delta_t * p.getV() + (delta_t * delta_t / (2 * p.getM())) * (p.getF());
    p.setX(new_x);
  });

  particle_container.single_iterator([this](Particle& p) {
    p.setOldF(p.getF());
    p.setF({0, 0, 0});
  });

  particle_container.pairIterator([this](Particle& p1, Particle& p2) {
    const dvec3 g12 = force.directionalForce(p1, p2);
    p1.setF(p1.getF() + g12); // this is basically F_i = \sum_j F_ij
    p2.setF(p2.getF() - g12); // this is because g12 = -g21
    // maybe we should only do one direction since other is covered before optimizing
    // also later on maybe force isn't bidirectional anymore? what then?
  });

  particle_container.single_iterator([this](Particle& p) {
    const dvec3 new_v = p.getV() + (delta_t / (2 * p.getM()) * (p.getOldF() + p.getF()));
    p.setV(new_v);
  });
}
